using FileIO
using Images
using Interpolations
using JLD2
using LinearAlgebra
using Plots
using Printf
using Random
using Statistics

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(PROJECT_ROOT, "src", "forwardmodel.jl"))
using .Model

# ── configuration ────────────────────────────────────────────────────────────

const ZOOM_LEVEL = "zoom5"
const INTRINSIC_DIR = joinpath(
    PROJECT_ROOT, "data", "calibration_camera", ZOOM_LEVEL, "detected",
)
const INTRINSIC_FILE = joinpath(INTRINSIC_DIR, "camint_all.jld2")
const INTRINSIC_COVARIANCE_FILE = joinpath(
    INTRINSIC_DIR, "subset_summaries", "reprojection_bootstrap_Nboot500.jld2",
)
const EXTRINSIC_DIR = joinpath(PROJECT_ROOT, "data", "test_droplets")
const EXTRINSIC_UNCERTAINTY_FILE = joinpath(
    EXTRINSIC_DIR, "extrinsic_uncertainty_$(ZOOM_LEVEL).jld2",
)
const ARUCO_IMAGE_FILE = joinpath(EXTRINSIC_DIR, "$(ZOOM_LEVEL).JPG")
const CORNER_RESULTS_FILE = joinpath(
    EXTRINSIC_DIR, "aruco_corner_results_$(ZOOM_LEVEL)_camtodrop.jld2",
)
const OUTPUT_DIR = joinpath(EXTRINSIC_DIR, "camtodrop_uncertainty")
const OUTPUT_FILE = joinpath(OUTPUT_DIR, "uncertainty_camtodrop.jld2")

const SIGMA_PIX = 0.3                 # px; isotropic contour-pixel noise
const SIGMA_POSE_OFFSET = 0.0         # mm; translation along the model-plane normal
const N_MODEL_POINTS = 5_000
const N_VIRTUAL_POINTS = 100          # Full covariance has size (2N)×(2N)
const N_LINEAR_MC = 500
const N_JOINT_MC = 2_000
const N_BASELINE_EXT_MC = 1_000
const MARKER_LENGTH_MM = 10.0
const SEED = 123

# Nominal pendant-drop and relative-pose values follow main_BI.jl.
const GAMMA_MN_PER_M = 72.0
const R_APEX_MM = 1.62
const R_NEEDLE_MM = 1.83 / 2
const THETA_X_DEG = 0.0
const THETA_Y_DEG = 0.0
const THETA_Z_DEG = 0.0
const T_REL_MM = Float64[-8.3, 0.6, 0.0]
const PLANE_COORDINATE_MM = 0.0

const CAMERA = Model.Camera
const UQ = CAMERA.IntrinsicUncertainty

function require_keys(data, keys, path)
    missing = filter(key -> !haskey(data, key), keys)
    isempty(missing) ||
        error("$path is missing required keys: $(join(missing, ", ")).")
end

function resample_curve(values::AbstractVector{<:Real}, count::Int)
    length(values) >= 2 || error("Physical contour has fewer than two points.")
    count >= 2 || throw(ArgumentError("Contour sample count must be at least two."))
    locations = range(1.0, Float64(length(values)); length=count)
    output = Vector{Float64}(undef, count)
    for (i, location) in enumerate(locations)
        left = floor(Int, location)
        right = min(left + 1, length(values))
        fraction = location - left
        output[i] = (1 - fraction) * values[left] + fraction * values[right]
    end
    return output
end

function stack_rz(r::AbstractVector, z::AbstractVector)
    length(r) == length(z) || throw(ArgumentError("r and z lengths differ."))
    output = Vector{promote_type(eltype(r), eltype(z))}(undef, 2 * length(r))
    for i in eachindex(r, z)
        output[2i - 1] = r[i]
        output[2i] = z[i]
    end
    return output
end

function unstack_rz(values::AbstractVector)
    iseven(length(values)) || throw(ArgumentError("Expected interleaved r,z values."))
    return values[1:2:end], values[2:2:end]
end

function contour_normals(r::AbstractVector, z::AbstractVector)
    length(r) == length(z) || throw(ArgumentError("r and z lengths differ."))
    n = length(r)
    n >= 2 || throw(ArgumentError("At least two contour points are required."))
    normals = Matrix{Float64}(undef, n, 2)
    for i in 1:n
        left, right = i == 1 ? (1, 2) : i == n ? (n - 1, n) : (i - 1, i + 1)
        dr = Float64(r[right] - r[left])
        dz = Float64(z[right] - z[left])
        magnitude = hypot(dr, dz)
        magnitude > 0 || error("Degenerate contour tangent at point $i.")
        normals[i, :] .= (-dz / magnitude, dr / magnitude)
    end
    return normals
end

function normal_covariance(covariance::AbstractMatrix, normals::AbstractMatrix)
    n = size(normals, 1)
    size(covariance) == (2n, 2n) ||
        throw(ArgumentError("Covariance does not match the contour length."))
    projected = Matrix{Float64}(undef, n, n)
    for i in 1:n, j in 1:n
        block = covariance[2i-1:2i, 2j-1:2j]
        projected[i, j] = dot(view(normals, i, :), block * view(normals, j, :))
    end
    return projected
end

normal_sigmas(covariance, normals) =
    sqrt.(max.(diag(normal_covariance(covariance, normals)), 0.0))

function block_diagonal_extrinsic(Σ_rvec, Σ_t)
    Σ = zeros(Float64, 6, 6)
    Σ[1:3, 1:3] .= Σ_rvec
    Σ[4:6, 4:6] .= Σ_t
    return Σ
end

function rotation_log(R::AbstractMatrix{<:Real})
    cosine = clamp((tr(R) - 1) / 2, -1.0, 1.0)
    θ = acos(cosine)
    vee = Float64[R[3, 2] - R[2, 3], R[1, 3] - R[3, 1], R[2, 1] - R[1, 2]]
    θ < 1e-8 && return vee ./ 2
    abs(sin(θ)) > 1e-12 ||
        error("Cannot compute a stable rotation-vector increment near angle π.")
    return (θ / (2sin(θ))) .* vee
end

function summarize_pose_samples(R_samples, t_samples)
    length(R_samples) == length(t_samples) ||
        throw(ArgumentError("Rotation and translation sample counts differ."))
    length(R_samples) >= 2 ||
        throw(ArgumentError("At least two pose samples are required."))

    R_mean = copy(first(R_samples))
    for _ in 1:10
        increments = [
            rotation_log(transpose(R_mean) * R) for R in R_samples
        ]
        mean_increment = vec(mean(reduce(hcat, increments); dims=2))
        R_mean = UQ._right_perturbed_rotation(R_mean, mean_increment)
        norm(mean_increment) < 1e-12 && break
    end

    samples = Matrix{Float64}(undef, length(R_samples), 6)
    for i in eachindex(R_samples, t_samples)
        samples[i, 1:3] .=
            rotation_log(transpose(R_mean) * R_samples[i])
        samples[i, 4:6] .= Float64.(t_samples[i])
    end
    t_mean = vec(mean(samples[:, 4:6]; dims=1))
    samples[:, 4:6] .-= reshape(t_mean, 1, 3)
    covariance = Matrix{Float64}(cov(samples; dims=1))
    return R_mean, t_mean, covariance
end

function correlation_matrix(covariance)
    standard_deviations = sqrt.(max.(diag(covariance), 0.0))
    denominator = standard_deviations * transpose(standard_deviations)
    return covariance ./ max.(denominator, eps(Float64))
end

function sample_covariance_from_pose_cache(R_mean, t_mean)
    candidate_files = [
        joinpath(EXTRINSIC_DIR, "extrinsic_bootstrap_$(ZOOM_LEVEL)_N1000.jld2"),
        joinpath(EXTRINSIC_DIR, "extrinsic_bootstrap_$(ZOOM_LEVEL)_N5000.jld2"),
        joinpath(
            EXTRINSIC_DIR,
            "extrinsic_bootstrap_$(ZOOM_LEVEL)_fit_corner_N1000.jld2",
        ),
    ]
    for path in candidate_files
        isfile(path) || continue
        data = JLD2.load(path)
        all(key -> haskey(data, key), ("R_samples", "t_samples")) || continue
        R_raw, t_raw = data["R_samples"], data["t_samples"]
        R_samples = R_raw isa AbstractArray{<:Real,3} ?
            [Matrix{Float64}(R_raw[i, :, :]) for i in axes(R_raw, 1)] :
            [Matrix{Float64}(R) for R in R_raw]
        t_samples = t_raw isa AbstractMatrix{<:Real} ?
            [Float64.(vec(t_raw[i, :])) for i in axes(t_raw, 1)] :
            [Float64.(t) for t in t_raw]
        length(R_samples) == length(t_samples) ||
            error("Mismatched pose samples in $path.")
        sample_matrix = Matrix{Float64}(undef, length(R_samples), 6)
        for i in eachindex(R_samples, t_samples)
            sample_matrix[i, 1:3] .= rotation_log(transpose(R_mean) * R_samples[i])
            sample_matrix[i, 4:6] .= t_samples[i] .- t_mean
        end
        @info "Using saved extrinsic bootstrap samples for full 6×6 covariance" path
        return Matrix{Float64}(cov(sample_matrix; dims=1)), path
    end
    @warn "No saved extrinsic pose samples found; using block-diagonal extrinsic covariance."
    return nothing, nothing
end

function load_or_fit_aruco_corners(K, dist)
    cached_paths = (
        joinpath(EXTRINSIC_DIR, "aruco_corner_results_$(ZOOM_LEVEL).jld2"),
        CORNER_RESULTS_FILE,
        joinpath(
            EXTRINSIC_DIR,
            "extrinsic_bootstrap_$(ZOOM_LEVEL)_fit_corner_N1000.jld2",
        ),
    )
    for path in cached_paths
        isfile(path) || continue
        cached = JLD2.load(path)
        haskey(cached, "Σ_corners") || continue
        corners_key = haskey(cached, "fitted_corners_px") ? "fitted_corners_px" :
                      haskey(cached, "corners_px") ? "corners_px" : nothing
        corners_key === nothing && continue
        haskey(cached, "marker_ids") || path != CORNER_RESULTS_FILE ||
            error("$path is missing marker_ids.")
        @info "Loading cached ArUco corners and covariance" path
        return Matrix{Float64}(cached[corners_key]),
               Matrix{Float64}[Matrix{Float64}(Σ) for Σ in cached["Σ_corners"]]
    end

    isfile(ARUCO_IMAGE_FILE) || error("ArUco image not found: $ARUCO_IMAGE_FILE")
    corners_px, marker_ids = UQ.detect_aruco_corners(ARUCO_IMAGE_FILE, K, dist)
    size(corners_px) == (4, 2) ||
        error("Expected four ArUco corners, received $(size(corners_px)).")

    image = FileIO.load(ARUCO_IMAGE_FILE)
    intensity = 1 .- Float64.(Gray.(image))
    interpolant = linear_interpolation(
        (1:size(intensity, 1), 1:size(intensity, 2)),
        intensity;
        extrapolation_bc=NaN,
    )
    fitted_corners = copy(corners_px)
    corner_results = Vector{Any}(undef, 4)
    for i in 1:4
        corner = corners_px[i, :] .+ 1.0
        previous = corners_px[mod1(i - 1, 4), :] .+ 1.0
        following = corners_px[mod1(i + 1, 4), :] .+ 1.0
        edge_a = previous .- corner
        edge_b = following .- corner
        length_a, length_b = norm(edge_a), norm(edge_b)
        shortest = min(length_a, length_b)
        max_distance = 0.4 * shortest
        max_distance > 13.0 ||
            error("ArUco corner $i cannot support the configured edge profiles.")
        result = UQ.fit_corner(
            interpolant,
            corner,
            edge_a ./ length_a,
            edge_b ./ length_b;
            profile_positions=collect(range(8.0, max_distance; length=6)),
            half=min(40.0, shortest / 12),
        )
        result === nothing && error("Could not fit ArUco corner $i.")
        corner_results[i] = result
        fitted_corners[i, :] .+= result.q
    end
    Σ_corners = Matrix{Float64}[Matrix{Float64}(result.Σ) for result in corner_results]
    mkpath(dirname(CORNER_RESULTS_FILE))
    JLD2.jldsave(
        CORNER_RESULTS_FILE;
        corner_results,
        opencv_corners_px=corners_px,
        fitted_corners_px=fitted_corners,
        Σ_corners,
        marker_ids,
        zoom_level=ZOOM_LEVEL,
    )
    @info "Saved regenerated ArUco observations and corner covariance" CORNER_RESULTS_FILE
    return fitted_corners, Σ_corners
end

function make_normal_band(r, z, normals, σ)
    upper = [Float64[r[i] + σ[i] * normals[i, 1],
                     z[i] + σ[i] * normals[i, 2]] for i in eachindex(r)]
    lower = [Float64[r[i] - σ[i] * normals[i, 1],
                     z[i] - σ[i] * normals[i, 2]] for i in eachindex(r)]
    polygon = vcat(upper, reverse(lower))
    return Shape(first.(polygon), last.(polygon))
end

function apex_radius_measure(r, z)
    indices = sortperm(abs.(r))[1:min(16, length(r))]
    length(indices) >= 4 || return NaN
    design = hcat(ones(length(indices)), Float64.(r[indices]).^2)
    coefficient = design \ Float64.(z[indices])
    curvature = abs(coefficient[2])
    return curvature > eps(Float64) ? 1 / (2curvature) : Inf
end

function covariance_terms(jacobians, Σ_intrinsics, Σ_ext_block, Σ_ext_full,
                          sigma_pix, sigma_pose_offset)
    Ji, Je, Jp, Jo = jacobians.J_intrinsics, jacobians.J_extrinsics,
                     jacobians.J_pixels, jacobians.J_offset
    camera_block = Ji * Σ_intrinsics * transpose(Ji) +
                   Je * Σ_ext_block * transpose(Je)
    camera_full = Ji * Σ_intrinsics * transpose(Ji) +
                  Je * Σ_ext_full * transpose(Je)
    pixel_term = (sigma_pix^2) .* (Jp * transpose(Jp))
    offset_term = (sigma_pose_offset^2) .* (Jo * transpose(Jo))
    return (
        camera_block=Matrix(Symmetric(camera_block)),
        camera_full=Matrix(Symmetric(camera_full)),
        pixel=Matrix(Symmetric(pixel_term)),
        offset=Matrix(Symmetric(offset_term)),
        total_block=Matrix(Symmetric(camera_block + pixel_term + offset_term)),
        total_full=Matrix(Symmetric(camera_full + pixel_term + offset_term)),
        camera_plus_offset=Matrix(Symmetric(camera_full + offset_term)),
    )
end

function sample_inverse_contour(
    pixels_u, pixels_v, phi, extrinsic_delta, pose_offset, R_mean, t_mean,
    R_rel, t_rel,
)
    parameters = vcat(phi, extrinsic_delta, pose_offset)
    r, z = UQ.inverse_project_camtodrop(
        pixels_u, pixels_v, parameters, R_mean, t_mean, R_rel, t_rel;
        plane_coordinate=PLANE_COORDINATE_MM,
    )
    return stack_rz(r, z)
end

function print_leading_modes(covariance, r, z; count=5)
    decomposition = eigen(Symmetric(covariance))
    order = sortperm(decomposition.values; rev=true)
    centroid = [mean(r), mean(z)]
    rigid = zeros(Float64, 2 * length(r), 3)
    for i in eachindex(r)
        rigid[2i - 1, 1] = 1
        rigid[2i, 2] = 1
        rigid[2i - 1, 3] = -(z[i] - centroid[2])
        rigid[2i, 3] = r[i] - centroid[1]
    end
    Qrigid = Matrix(qr(rigid).Q)[:, 1:3]
    scale = stack_rz(r .- centroid[1], z .- centroid[2])
    scale_residual = scale - Qrigid * (transpose(Qrigid) * scale)
    scale_residual ./= max(norm(scale_residual), eps(Float64))
    println("\nLeading camera-covariance modes (translation/rotation, scale, shape):")
    rows = min(count, length(order))
    for rank in 1:rows
        index = order[rank]
        vector = decomposition.vectors[:, index]
        rigid_fraction = sum(abs2, transpose(Qrigid) * vector)
        scale_fraction = abs2(dot(scale_residual, vector))
        shape_fraction = max(0.0, 1 - rigid_fraction - scale_fraction)
        @printf(
            "  %d: λ=%.6g  rigid=%.1f%%  scale=%.1f%%  shape=%.1f%%\n",
            rank, decomposition.values[index], 100 * rigid_fraction,
            100 * scale_fraction, 100 * shape_fraction,
        )
    end
    return decomposition.values[order], decomposition.vectors[:, order]
end

function main()
    SIGMA_PIX >= 0 || error("SIGMA_PIX must be nonnegative.")
    SIGMA_POSE_OFFSET >= 0 || error("SIGMA_POSE_OFFSET must be nonnegative.")
    N_JOINT_MC >= 2 && N_LINEAR_MC >= 2 ||
        error("Monte Carlo sample counts must be at least two.")

    for path in (
        INTRINSIC_FILE, INTRINSIC_COVARIANCE_FILE,
        EXTRINSIC_UNCERTAINTY_FILE,
    )
        isfile(path) || error("Required calibration file not found: $path")
    end
    intrinsics_file = JLD2.load(INTRINSIC_FILE)
    covariance_file = JLD2.load(INTRINSIC_COVARIANCE_FILE)
    extrinsics_file = JLD2.load(EXTRINSIC_UNCERTAINTY_FILE)
    require_keys(intrinsics_file, ("camera_matrix", "dist_coeffs"), INTRINSIC_FILE)
    require_keys(covariance_file, ("Σ_φ",), INTRINSIC_COVARIANCE_FILE)
    require_keys(
        extrinsics_file, ("rvec_mean", "t_mean", "Σ_rvec", "Σ_t"),
        EXTRINSIC_UNCERTAINTY_FILE,
    )

    K = Matrix{Float64}(intrinsics_file["camera_matrix"])
    dist = vec(Float64.(intrinsics_file["dist_coeffs"]))[1:5]
    intrinsics = Float64[
        K[1, 1], K[2, 2], K[1, 3], K[2, 3], dist[1], dist[2],
        dist[3], dist[4], dist[5],
    ]
    Σ_intrinsics = Matrix{Float64}(covariance_file["Σ_φ"])
    size(Σ_intrinsics) == (9, 9) ||
        error("Σ_φ must be 9×9 in [fx,fy,cx,cy,k1,k2,p1,p2,k3] order.")
    legacy_Σ_rvec = Matrix{Float64}(extrinsics_file["Σ_rvec"])
    legacy_Σ_t = Matrix{Float64}(extrinsics_file["Σ_t"])
    fitted_corners, Σ_corners = load_or_fit_aruco_corners(K, dist)
    baseline_intrinsic_samples =
        repeat(reshape(intrinsics, 1, 9), N_BASELINE_EXT_MC, 1)
    baseline_R_samples, baseline_t_samples =
        UQ.joint_intrinsic_extrinsic_bootstrap(
            fitted_corners, Σ_corners, baseline_intrinsic_samples,
            MARKER_LENGTH_MM;
            seed=SEED + 3,
        )
    R_mean, t_mean, Σ_ext_full =
        summarize_pose_samples(baseline_R_samples, baseline_t_samples)
    Σ_rvec = Matrix(Σ_ext_full[1:3, 1:3])
    Σ_t = Matrix(Σ_ext_full[4:6, 4:6])
    Σ_ext_block = block_diagonal_extrinsic(Σ_rvec, Σ_t)
    extrinsic_sample_file = "recomputed from fitted corners and their covariance"

    legacy_Σ_ext = block_diagonal_extrinsic(legacy_Σ_rvec, legacy_Σ_t)
    stored_sd = sqrt.(max.(diag(legacy_Σ_ext), 0.0))
    baseline_sd = sqrt.(max.(diag(Σ_ext_full), 0.0))
    println("\nFixed-intrinsics bootstrap sanity check (recomputed / stored SD):")
    println("  ratios [δrvec(rad), δt(mm)] = ",
            baseline_sd ./ max.(stored_sd, eps(Float64)))
    println("  recomputed 6×6 pose correlation matrix:")
    display(round.(correlation_matrix(Σ_ext_full); digits=3))
    println("  stored rotation correlation matrix:")
    display(round.(correlation_matrix(legacy_Σ_rvec); digits=3))
    println("  stored translation correlation matrix:")
    display(round.(correlation_matrix(legacy_Σ_t); digits=3))
    # Bootstrap increments are right perturbations: R_sample = R_mean * exp(δrvec).

    physical_r, physical_z, _, hits = Model.PD_fn(
        log.([GAMMA_MN_PER_M, R_APEX_MM]),
        R_NEEDLE_MM;
        n_points=N_MODEL_POINTS,
        r_tol=1e-12,
        h_tol=1e-12,
        hits_max=2,
    )
    isempty(physical_r) && error("PD_fn returned an empty contour.")
    hits > 0 || @warn "PD_fn detected no needle-radius crossings" hits
    contour_r = resample_curve(physical_r, N_VIRTUAL_POINTS)
    contour_z = resample_curve(physical_z, N_VIRTUAL_POINTS)
    R_rel = CAMERA.getRrel_fn(THETA_X_DEG, THETA_Y_DEG, THETA_Z_DEG)
    base_parameters = vcat(intrinsics, zeros(6), 0.0)

    virtual_u, virtual_v, _ = UQ.forward_project_camtodrop(
        contour_r, contour_z, base_parameters,
        R_mean, t_mean, R_rel, T_REL_MM;
        plane_coordinate=PLANE_COORDINATE_MM,
    )
    recovered_r, recovered_z = UQ.inverse_project_camtodrop(
        virtual_u, virtual_v, base_parameters,
        R_mean, t_mean, R_rel, T_REL_MM;
        plane_coordinate=PLANE_COORDINATE_MM,
    )
    roundtrip_errors = hypot.(recovered_r .- contour_r, recovered_z .- contour_z)
    max_roundtrip_error = maximum(roundtrip_errors)
    println("Maximum forward/inverse round-trip error: $(max_roundtrip_error) mm")

    jacobians = UQ.camtodrop_jacobians(
        virtual_u, virtual_v, intrinsics,
        R_mean, t_mean, R_rel, T_REL_MM;
        plane_coordinate=PLANE_COORDINATE_MM,
    )
    terms = covariance_terms(
        jacobians, Σ_intrinsics, Σ_ext_block, Σ_ext_full,
        SIGMA_PIX, SIGMA_POSE_OFFSET,
    )
    normals = contour_normals(jacobians.coordinates[1:2:end],
                              jacobians.coordinates[2:2:end])
    independent_normal_sigma = normal_sigmas(terms.camera_full + terms.offset, normals)

    rng = MersenneTwister(SEED)
    linear_samples = Matrix{Float64}(undef, N_LINEAR_MC, 2N_VIRTUAL_POINTS)
    Σ_independent = zeros(Float64, 15, 15)
    Σ_independent[1:9, 1:9] .= Σ_intrinsics
    Σ_independent[10:15, 10:15] .= Σ_ext_full
    parameter_mean = vcat(intrinsics, zeros(6))
    for sample_index in 1:N_LINEAR_MC
        sampled_parameters = UQ.sample_gaussian(rng, parameter_mean, Σ_independent)
        offset_sample = SIGMA_POSE_OFFSET * randn(rng)
        linear_samples[sample_index, :] .= sample_inverse_contour(
            virtual_u, virtual_v, sampled_parameters[1:9],
            sampled_parameters[10:15], offset_sample,
            R_mean, t_mean, R_rel, T_REL_MM,
        )
    end
    Σ_linear_mc = Matrix{Float64}(cov(linear_samples; dims=1))
    mc_normal_sigma = normal_sigmas(Σ_linear_mc, normals)
    linear_normal_sigma = normal_sigmas(terms.camera_full + terms.offset, normals)
    println("\nCovariance component summaries [model-local mm²]:")
    for (label, matrix) in (
        ("camera, block-diagonal extrinsics", terms.camera_block),
        ("camera, full 6×6 extrinsics", terms.camera_full),
        ("pixel noise", terms.pixel),
        ("relative-pose normal offset", terms.offset),
        ("total, full 6×6 extrinsics", terms.total_full),
    )
        println("  $label: trace=$(tr(matrix))")
    end
    println(
        "  median normal σ (linear / MC): ",
        "$(median(linear_normal_sigma)) / $(median(mc_normal_sigma)) mm",
    )

    values_linear, vectors_linear =
        print_leading_modes(terms.camera_full, jacobians.coordinates[1:2:end],
                            jacobians.coordinates[2:2:end])
    values_mc, vectors_mc =
        print_leading_modes(Σ_linear_mc, jacobians.coordinates[1:2:end],
                            jacobians.coordinates[2:2:end])
    println("\nLeading covariance eigenvalue comparison (linear vs MC):")
    for i in 1:min(5, length(values_linear), length(values_mc))
        cosine = abs(dot(vectors_linear[:, i], vectors_mc[:, i]))
        @printf("  %d: %.6g / %.6g, |dot(vectors)|=%.4f\n",
                i, values_linear[i], values_mc[i], cosine)
    end

    intrinsic_samples = Matrix{Float64}(undef, N_JOINT_MC, 9)
    rng_intrinsic = MersenneTwister(SEED + 1)
    for i in 1:N_JOINT_MC
        intrinsic_samples[i, :] .=
            UQ.sample_gaussian(rng_intrinsic, intrinsics, Σ_intrinsics)
    end
    R_joint, t_joint = UQ.joint_intrinsic_extrinsic_bootstrap(
        fitted_corners, Σ_corners, intrinsic_samples, MARKER_LENGTH_MM;
        seed=SEED + 2,
    )
    joint_camera_samples = Matrix{Float64}(undef, N_JOINT_MC, 2 * N_VIRTUAL_POINTS)
    joint_samples = Matrix{Float64}(undef, N_JOINT_MC, 2N_VIRTUAL_POINTS)
    joint_apex_radius = Vector{Float64}(undef, N_JOINT_MC)
    for i in 1:N_JOINT_MC
        rotation_delta = rotation_log(transpose(R_mean) * R_joint[i])
        translation_delta = t_joint[i] .- t_mean
        offset_sample = SIGMA_POSE_OFFSET * randn(rng_intrinsic)
        joint_camera_samples[i, :] .= sample_inverse_contour(
            virtual_u, virtual_v, intrinsic_samples[i, :],
            vcat(rotation_delta, translation_delta), 0.0,
            R_mean, t_mean, R_rel, T_REL_MM,
        )
        joint_samples[i, :] .= sample_inverse_contour(
            virtual_u, virtual_v, intrinsic_samples[i, :],
            vcat(rotation_delta, translation_delta), offset_sample,
            R_mean, t_mean, R_rel, T_REL_MM,
        )
        r_sample, z_sample = unstack_rz(view(joint_samples, i, :))
        joint_apex_radius[i] = apex_radius_measure(r_sample, z_sample)
    end
    Σ_joint_camera_mc = Matrix{Float64}(cov(joint_camera_samples; dims=1))
    Σ_joint_mc = Matrix{Float64}(cov(joint_samples; dims=1))
    Σ_total_joint = Matrix(Symmetric(Σ_joint_mc + terms.pixel))
    joint_normal_sigma = normal_sigmas(Σ_joint_mc, normals)
    joint_camera_normal_sigma = normal_sigmas(Σ_joint_camera_mc, normals)
    joint_total_normal_sigma = normal_sigmas(Σ_total_joint, normals)
    pixel_normal_sigma = normal_sigmas(terms.pixel, normals)
    normal_ratio = joint_normal_sigma ./ max.(independent_normal_sigma, eps(Float64))
    println("\nJoint-bootstrap comparison (joint reference / independent model):")
    println(
        "  median pointwise normal-σ ratio: $(round(median(normal_ratio), digits=3))",
    )
    println(
        "  10–90% normal-σ ratio: ",
        "$(round(quantile(normal_ratio, 0.1), digits=3))–",
        "$(round(quantile(normal_ratio, 0.9), digits=3))",
    )

    independent_apex_radius =
        [apex_radius_measure(view(linear_samples, i, 1:2:size(linear_samples, 2)),
                             view(linear_samples, i, 2:2:size(linear_samples, 2)))
         for i in 1:N_LINEAR_MC]
    size_ratio = std(joint_apex_radius) / max(std(independent_apex_radius), eps(Float64))
    println(
        "  apex-radius spread, independent / joint: ",
        "$(round(std(independent_apex_radius), digits=5)) / ",
        "$(round(std(joint_apex_radius), digits=5)) mm",
    )
    println("  joint / independent apex-radius σ ratio: $(round(size_ratio, digits=3))")
    median(normal_ratio) < 0.8 ?
        println("  Independence substantially inflates normal uncertainty for this contour.") :
        median(normal_ratio) > 1.2 ?
        println("  Independence materially underestimates normal uncertainty for this contour.") :
        println("  No substantial median normal-uncertainty inflation is indicated.")

    mkpath(OUTPUT_DIR)
    plot_linear_mc = plot(
        1:N_VIRTUAL_POINTS, linear_normal_sigma;
        label="Linearized", linewidth=2, xlabel="Contour point",
        ylabel="Normal-direction σ [mm]",
        title="Linearized vs Monte Carlo uncertainty",
    )
    plot!(plot_linear_mc, 1:N_VIRTUAL_POINTS, mc_normal_sigma;
          label="Independent Monte Carlo", linewidth=2, linestyle=:dash)
    plot!(plot_linear_mc, 1:N_VIRTUAL_POINTS, joint_normal_sigma;
          label="Joint bootstrap camera", linewidth=2, linestyle=:dot)
    plot!(plot_linear_mc, 1:N_VIRTUAL_POINTS, joint_total_normal_sigma;
          label="Joint bootstrap total (camera + pixel)", linewidth=2, linestyle=:solid)
    plot!(plot_linear_mc, 1:N_VIRTUAL_POINTS, joint_camera_normal_sigma;
          label="Joint camera term", linewidth=2, linestyle=:dashdot)
    plot!(plot_linear_mc, 1:N_VIRTUAL_POINTS, pixel_normal_sigma;
          label="Pixel noise (σ = $(SIGMA_PIX) px)", linewidth=2, linestyle=:dash)
    savefig(plot_linear_mc, joinpath(OUTPUT_DIR, "linear_vs_monte_carlo.png"))

    plot_bands = plot(
        jacobians.coordinates[1:2:end], jacobians.coordinates[2:2:end];
        color=:black, linewidth=2, label="Virtual contour",
        xlabel="r [mm]", ylabel="z [mm]",
        title="Normal-direction uncertainty bands",
        aspect_ratio=:equal,
    )
    plot!(plot_bands, make_normal_band(
        jacobians.coordinates[1:2:end], jacobians.coordinates[2:2:end],
        normals, independent_normal_sigma,
    ); color=:blue, fillalpha=0.16, linealpha=0, label="Independent, 1σ")
    plot!(plot_bands, make_normal_band(
        jacobians.coordinates[1:2:end], jacobians.coordinates[2:2:end],
        normals, joint_normal_sigma,
    ); color=:orange, fillalpha=0.16, linealpha=0, label="Joint bootstrap, 1σ")
    savefig(plot_bands, joinpath(OUTPUT_DIR, "normal_uncertainty_bands.png"))

    Σ_normal_camera = normal_covariance(terms.camera_full, normals)
    Σ_normal_pixel = normal_covariance(terms.pixel, normals)
    Σ_normal_total = normal_covariance(terms.total_full, normals)
    Σ_normal_joint = normal_covariance(Σ_joint_mc, normals)
    metadata = Dict(
        "parameter_order" => collect(jacobians.parameter_order),
        "intrinsic_order" => [
            "fx_px", "fy_px", "cx_px", "cy_px", "k1", "k2", "p1", "p2", "k3",
        ],
        "extrinsic_order" => [
            "δrx_right_rad", "δry_right_rad", "δrz_right_rad",
            "δtx_mm", "δty_mm", "δtz_mm",
        ],
        "units" => Dict("contour" => "mm", "pixel" => "px",
                        "rotation" => "rad", "covariance" => "mm^2"),
        "frame_definition" =>
            "Model-local meridian (r,z), third coordinate plane_coordinate=0; mm.",
        "extrinsic_rotation_perturbation" =>
            "Right perturbation: R_sample = R_mean * exp(δrvec); δrvec = log(R_mean' * R_sample).",
        "independence_assumption" =>
            "Intrinsic calibration and extrinsic pose uncertainty are independent. Both independent and joint calculations use extrinsic uncertainty generated from the same fitted corner locations and fitted corner covariances.",
        "extrinsic_noise_basis" =>
            "Each fitted corner is perturbed independently using its fitted 2×2 covariance; nominal covariance and joint bootstrap use the same centers and covariance matrices.",
        "virtual_plane_assumption" =>
            "The PD_fn contour is forward-projected with its third model coordinate fixed at zero so it lies on the inverse-projection plane.",
        "extrinsic_sample_source" => extrinsic_sample_file,
        "corner_results_file" => CORNER_RESULTS_FILE,
        "baseline_extrinsic_bootstrap_samples" => N_BASELINE_EXT_MC,
        "joint_bootstrap_samples" => N_JOINT_MC,
        "baseline_to_stored_extrinsic_sd_ratio" =>
            baseline_sd ./ max.(stored_sd, eps(Float64)),
        "seed" => SEED,
        "sigma_pix_px" => SIGMA_PIX,
        "sigma_pose_offset_mm" => SIGMA_POSE_OFFSET,
        "roundtrip_max_error_mm" => max_roundtrip_error,
    )
    JLD2.jldsave(
        OUTPUT_FILE;
        virtual_r=jacobians.coordinates[1:2:end],
        virtual_z=jacobians.coordinates[2:2:end],
        virtual_u=Float64.(virtual_u),
        virtual_v=Float64.(virtual_v),
        roundtrip_errors,
        Σ_intrinsics,
        Σ_rvec,
        Σ_t,
        rvec_mean=rotation_log(R_mean),
        t_mean,
        Σ_ext_block,
        Σ_ext_full,
        Σ_ext_legacy_blockdiag=legacy_Σ_ext,
        baseline_to_stored_extrinsic_sd_ratio=
            baseline_sd ./ max.(stored_sd, eps(Float64)),
        Σ_camera_independent_block=terms.camera_block,
        Σ_camera_independent_full_ext=terms.camera_full,
        Σ_pixel=terms.pixel,
        Σ_pose_offset=terms.offset,
        Σ_total_independent_block=terms.total_block,
        Σ_total_independent_full_ext=terms.total_full,
        Σ_linear_monte_carlo=Σ_linear_mc,
        Σ_camera_joint_reference=Σ_joint_camera_mc,
        Σ_pixel_joint_reference=terms.pixel,
        Σ_total_joint_reference=Σ_total_joint,
        Σ_joint_monte_carlo=Σ_joint_mc,
        Σ_normal_camera=Σ_normal_camera,
        Σ_normal_pixel=Σ_normal_pixel,
        Σ_normal_total=Σ_normal_total,
        Σ_normal_joint=Σ_normal_joint,
        normal_sigma_independent=independent_normal_sigma,
        normal_sigma_linearized=linear_normal_sigma,
        normal_sigma_monte_carlo=mc_normal_sigma,
        normal_sigma_joint=joint_normal_sigma,
        normal_sigma_joint_camera=joint_camera_normal_sigma,
        normal_sigma_joint_total=joint_total_normal_sigma,
        normal_sigma_pixel=pixel_normal_sigma,
        normal_sigma_ratio=normal_ratio,
        joint_apex_radius_mm=joint_apex_radius,
        independent_apex_radius_mm=independent_apex_radius,
        metadata,
        fitted_aruco_corners_px=fitted_corners,
        Σ_corners,
        baseline_R_samples,
        baseline_t_samples,
    )
    println("\nSaved virtual contour and covariance results: $OUTPUT_FILE")
    println("Saved plots in: $OUTPUT_DIR")
end

if abspath(PROGRAM_FILE) == (@__FILE__) || isinteractive()
    main()
end
