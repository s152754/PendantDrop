"""
Camera calibration utilities

This include-only file provides checkerboard image selection and calibration,
intrinsic-parameter uncertainty estimation, and plotting helpers. It relies on
packages and helpers loaded by the parent camera module; do not include it alone.
"""


# Table of contents
# 1. Helper functions
# 2. Intrinsic calibration uncertainty
# 3. ArUco corner uncertainty
# 4. Extrinsic pose uncertainty
# 5. Camera-to-drop projection
# 6. Uncertainty propagation and sampling
# 7. Shape and position analysis
# 8. Not uncertainty-related

using Printf
using Images: Gray
using Interpolations: linear_interpolation
import ForwardDiff
const Camera = parentmodule(@__MODULE__)

# =============================================================================
# Helper functions
# Shared conversions, validation, contour packing, and covariance utilities.
# =============================================================================
# Internal conversion helpers
"""
Convert a camera matrix to a 3×3 `Float64` matrix.

input:
- Kmat [matrix, nested vector, or `nothing`]

outputs:
- Camera matrix [unitless]
"""
function _coerce_camera_matrix(Kmat)
    Kmat === nothing && return zeros(3, 3)
    Kmat isa AbstractMatrix && return Float64.(Matrix(Kmat))
    if Kmat isa AbstractVector
        length(Kmat) == 9 && return reshape(Float64.(Kmat), 3, 3)
        if !isempty(Kmat) && all(x -> x isa AbstractVector || x isa AbstractMatrix, Kmat)
            rows = [reshape(Float64.(collect(row)), 1, length(row)) for row in Kmat]
            return reduce(vcat, rows)
        end
    end
    return Float64.(Matrix(Kmat))
end

"""
Convert distortion coefficients to a flat `Float64` vector.

input:
- distcoef [coefficient vector, matrix, scalar, or `nothing`]

outputs:
- Distortion coefficients [unitless]
"""
function _coerce_distortion_coeffs(distcoef)
    distcoef === nothing && return Float64[]
    distcoef isa AbstractMatrix && return vec(Float64.(distcoef))
    if distcoef isa AbstractVector
        !isempty(distcoef) && distcoef[1] isa AbstractVector &&
            return Float64.(collect(Iterators.flatten(distcoef)))
        return Float64.(distcoef)
    end
    return Float64[Float64(distcoef)]
end

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

function block_diagonal_extrinsic(Σ_rvec, Σ_t)
    Σ = zeros(Float64, 6, 6)
    Σ[1:3, 1:3] .= Σ_rvec
    Σ[4:6, 4:6] .= Σ_t
    return Σ
end

function correlation_matrix(covariance)
    standard_deviations = sqrt.(max.(diag(covariance), 0.0))
    denominator = standard_deviations * transpose(standard_deviations)
    return covariance ./ max.(denominator, eps(Float64))
end

function _validate_correlation_matrix(correlation, n::Int)
    size(correlation) == (n, n) ||
        throw(ArgumentError("Correlation matrix must be $n×$n."))
    all(isfinite, correlation) ||
        throw(ArgumentError("Correlation matrix must contain finite values."))
    isapprox(correlation, transpose(correlation); atol=1e-10, rtol=1e-8) ||
        throw(ArgumentError("Correlation matrix must be symmetric."))
    isapprox(diag(correlation), ones(n); atol=1e-10, rtol=0) ||
        throw(ArgumentError("Correlation matrix must have a unit diagonal."))
    symmetric = Matrix(Symmetric(Float64.(correlation)))
    tolerance = 1e-10 * max(opnorm(symmetric, Inf), 1.0)
    minimum(eigvals(Symmetric(symmetric))) >= -tolerance ||
        throw(ArgumentError("Correlation matrix must be positive semidefinite."))
    return symmetric
end

# =============================================================================
# Intrinsic calibration uncertainty
# Residual covariance, correlated calibration bootstrap, and subset variability.
# =============================================================================


# Reprojection and corner-noise uncertainty

"""
Compute per-corner reprojection residuals and their x/y standard deviations.

input:
- objpoints [3-D points]
- imgpoints [observed image points, px]
- rvecs [rotation vectors]
- tvecs [translation vectors]
- K [camera matrix]
- dist [distortion coefficients]

outputs:
- sigma_x, sigma_y [px]
- residuals [N×2 array, px]
- Sigma_corner [2×2 covariantiematrix corner-ruis, px²]
"""
function compute_reprojection_residuals(objpoints, imgpoints, rvecs, tvecs, K, dist)
    py"""
import cv2, numpy as np
all_residuals = []
for objp, imgp, rvec, tvec in zip($objpoints, $imgpoints, $rvecs, $tvecs):
    objp = np.array(objp, dtype=np.float32)
    imgp = np.array(imgp, dtype=np.float64)
    rvec = np.array(rvec, dtype=np.float64).reshape(3,1)
    tvec = np.array(tvec, dtype=np.float64).reshape(3,1)
    projected, _ = cv2.projectPoints(objp, rvec, tvec,
                                     np.array($K, dtype=np.float64),
                                     np.array($dist, dtype=np.float64))
    all_residuals.append(imgp - projected.reshape(-1,2))
all_residuals = np.vstack(all_residuals)
Sigma_corner = np.cov(all_residuals[:,0],
                      all_residuals[:,1]).tolist()
result = (float(np.std(all_residuals[:,0])),
          float(np.std(all_residuals[:,1])),
          all_residuals.tolist(),
          Sigma_corner)
"""
    sigma_x, sigma_y, residuals, Sigma_corner = py"result"
    return sigma_x, sigma_y, residuals, Sigma_corner
end


"""
Estimate intrinsic-parameter uncertainty by perturbing corners and recalibrating.

input:
- objpoints [3-D points]
- imgpoints [detected image corners, px]
- image_size [width, height, px]
- Sigma_corner [2×2 covariantiematrix corner-ruis, px²]
- N_boot [bootstrap runs]

outputs:
- samples [intrinsic parameter samples]
- mean_phi [mean of fx, fy, cx, cy, k1, k2, p1, p2, k3]
- cov_phi [parameter covariance]
"""
function bootstrap_calibration(objpoints, imgpoints,
                               image_size, Sigma_corner::Matrix{Float64};
                               N_boot::Int = 500)
    py"""
import cv2, numpy as np
objpoints  = $objpoints
imgpoints  = $imgpoints
image_size = tuple($image_size)
Sigma      = np.array($Sigma_corner)
N_boot     = int($N_boot)
samples    = []
for i in range(N_boot):
    perturbed = []
    for corners in imgpoints:
        corners = np.array(corners, dtype=np.float32).reshape(-1,2)
        n_corners = len(corners)
        noise = np.random.multivariate_normal(
            [0, 0], Sigma, n_corners)
        noise_x = noise[:, 0]
        noise_y = noise[:, 1]
        noise = np.column_stack([noise_x, noise_y])
        perturbed.append((corners + noise).astype(np.float32).reshape(-1,1,2))
    objpoints_cv = [np.array(o, dtype=np.float32) for o in objpoints]
    ret, K, dist, _, _ = cv2.calibrateCamera(objpoints_cv, perturbed, image_size, None, None)
    dist = np.array(dist).reshape(-1)
    samples.append([K[0,0], K[1,1], K[0,2], K[1,2]] + dist.tolist())
    if (i + 1) % 50 == 0 or (i + 1) == N_boot:
        print(f"\r  Bootstrap: {i+1}/{N_boot} samples", flush=True)
print()
samples  = np.array(samples)
result   = (samples.tolist(), np.mean(samples, axis=0).tolist(),
            np.cov(samples, rowvar=False).tolist())
"""
    return py"result"
end


# Uncertainty from calibration image count

"""
Run and cache repeated calibrations over random image subsets.

input:
- img_dir [directory path]
- out_dir [reference-calibration directory path]
- summary_dir [subset-cache directory path]
- subset_sizes [images per subset]
- n_seeds [runs per subset size]

outputs:
- results [Dict mapping subset size to calibration samples]
- K_ref [reference camera matrix]
- dist_ref [reference distortion coefficients]
"""
function run_subset_study(
    img_dir::String,
    out_dir::String,
    summary_dir::String;
    subset_sizes::Vector{Int} = [5, 10, 20, 30, 40, 45],
    n_seeds::Int              = 10,
)
    mkpath(summary_dir)

    # reference calibration (all images)
    all_ref_path = joinpath(out_dir, "camint_all.jld2")
    if isfile(all_ref_path)
        @info "Loading reference calibration from $all_ref_path"
        ref      = load(all_ref_path)
        K_ref    = _coerce_camera_matrix(
            get(ref, "camera_matrix", get(ref, "K", get(ref, "K_all", nothing))),
        )
        dist_ref = _coerce_distortion_coeffs(
            get(ref, "dist_coeffs", get(ref, "dist", get(ref, "distcoef", get(ref, "distcoef_all", nothing)))),
        )
    else
        @info "Running reference calibration on all images …"
        all_files = list_camera_images(img_dir)
        isempty(all_files) && error("No calibration images found in $img_dir")
        ref_result = calibrate_camera_from_images(all_files)
        get(ref_result, "success", false) || error(
            "Reference camera calibration failed for images in $img_dir; " *
            "check the checkerboard settings and image set.",
        )
        K_ref      = _coerce_camera_matrix(ref_result["camera_matrix"])
        dist_ref   = _coerce_distortion_coeffs(ref_result["dist_coeffs"])
        objpoints = ref_result["objpoints"]
        imgpoints = ref_result["imgpoints"]
        rvecs = ref_result["rvecs"]
        tvecs = ref_result["tvecs"]
        image_size = ref_result["image_size"]
        jldsave(all_ref_path; camera_matrix=K_ref, dist_coeffs=dist_ref,
            rms=ref_result["rms"], objpoints, imgpoints, rvecs, tvecs, image_size)
    end

    if size(K_ref) != (3, 3) || length(dist_ref) < 5 ||
       !all(isfinite, K_ref) || !all(isfinite, dist_ref[1:5])
        error("Reference calibration at $all_ref_path has an invalid camera matrix " *
              "or fewer than five finite distortion coefficients; remove or regenerate it.")
    end

    results = Dict{Int,Dict}()

    for n in subset_sizes
        cache = joinpath(summary_dir, "subset_$(n)_seeds$(n_seeds).jld2")
        if isfile(cache)
            @info "Loading cached subset n=$n from $cache"
            d = load(cache)
            seed_ids = get(d, "seed_ids", nothing)
            if seed_ids === nothing
                if length(d["K_samples"]) == n_seeds
                    seed_ids = collect(1:n_seeds)
                else
                    @warn "Cached subset n=$n predates seed ID storage; plotting successful-run indices instead."
                    seed_ids = collect(1:length(d["K_samples"]))
                end
            end
            results[n] = Dict("K_samples"    => d["K_samples"],
                              "dist_samples" => d["dist_samples"],
                              "rms_values"   => d["rms_values"],
                              "seed_ids"     => seed_ids)
            continue
        end

        @info "Subset study  n=$n  ($n_seeds seeds) …"
        K_samples    = Matrix{Float64}[]
        dist_samples = Vector{Float64}[]
        rms_values   = Float64[]
        seed_ids     = Int[]

        for seed in 1:n_seeds
            files = random_image_subset(img_dir; subset_size=n, seed=seed)
            isempty(files) && continue
            res = calibrate_camera_from_images(files)
            res["success"] || continue
            push!(K_samples,    _coerce_camera_matrix(res["camera_matrix"]))
            push!(dist_samples, _coerce_distortion_coeffs(res["dist_coeffs"]))
            push!(seed_ids, seed)
            isnothing(res["rms"]) || push!(rms_values, Float64(res["rms"]))
        end

        if isempty(K_samples)
            @warn "No successful calibration runs for n=$n"
            continue
        end

        results[n] = Dict("K_samples" => K_samples, "dist_samples" => dist_samples,
                  "rms_values" => rms_values, "seed_ids" => seed_ids)
        jldsave(cache; K_samples, dist_samples, rms_values, seed_ids)
        @info "  Saved → $cache"
    end

    return results, K_ref, dist_ref
end


"""
Estimate intrinsic-parameter uncertainty by bootstrapping reprojection corner noise.

input:
- out_dir [reference-calibration directory path]
- summary_dir [bootstrap-cache and plot directory path]
- N_boot [bootstrap runs, at least 2]

outputs:
- sigma_x, sigma_y [corner-noise standard deviations, px]
- Sigma_corner [2×2 covariantiematrix corner-ruis, px²]
- residuals [N×2 reprojection residuals, px; returned for plotting, not cached]
- samples [N_boot×9 intrinsic-parameter samples]
- mean_φ [mean intrinsic-parameter vector]
- Σ_φ [9×9 covariance matrix]
- N_boot [bootstrap runs]
"""
function run_reprojection_bootstrap(
    out_dir::String,
    summary_dir::String;
    N_boot::Int = 500,
)
    N_boot >= 2 || throw(ArgumentError("N_boot must be at least 2 to estimate covariance."))

    reference_path = joinpath(out_dir, "camint_all.jld2")
    isfile(reference_path) || error(
        "Reference calibration not found at $reference_path. Run run_subset_study first.",
    )
    reference = load(reference_path)
    required_keys = ("camera_matrix", "dist_coeffs", "rms", "objpoints", "imgpoints",
                     "rvecs", "tvecs", "image_size")
    missing_keys = filter(key -> !haskey(reference, key), required_keys)
    isempty(missing_keys) || error(
        "Reference calibration $reference_path is missing required keys: " *
        "$(join(missing_keys, ", ")). Delete camint_all.jld2 and re-run " *
        "run_subset_study to regenerate it with reprojection data.",
    )

    K = _coerce_camera_matrix(reference["camera_matrix"])
    dist = _coerce_distortion_coeffs(reference["dist_coeffs"])
    sigma_x, sigma_y, residuals, Sigma_corner =
        compute_reprojection_residuals(
        reference["objpoints"], reference["imgpoints"], reference["rvecs"],
        reference["tvecs"], K, dist,
    )
    plot_reprojection_residuals(residuals; save_dir=summary_dir)

    cache_path = joinpath(summary_dir, "reprojection_bootstrap_Nboot$(N_boot).jld2")
    if isfile(cache_path)
        @info "Loading cached reprojection bootstrap from $cache_path"
        cached = load(cache_path)
        return sigma_x, sigma_y, residuals, Sigma_corner, cached["samples"],
               cached["mean_φ"], cached["Σ_φ"], Int(cached["N_boot"])
    end

    samples_raw, _, _ = bootstrap_calibration(
        reference["objpoints"], reference["imgpoints"], reference["image_size"],
        Sigma_corner; N_boot=N_boot,
    )
    samples = Matrix{Float64}(samples_raw)
    size(samples) == (N_boot, 9) || error(
        "Bootstrap returned samples with size $(size(samples)); expected ($N_boot, 9).",
    )
    mean_φ = vec(mean(samples; dims=1))
    Σ_φ = cov(samples; dims=1)

    mkpath(summary_dir)
    jldsave(cache_path; sigma_x, sigma_y, Sigma_corner, residuals, samples, mean_φ, Σ_φ, N_boot)
    @info "Saved reprojection bootstrap to $cache_path"
    return sigma_x, sigma_y, residuals, Sigma_corner, samples, mean_φ, Σ_φ, N_boot
end


"""
Summarize the mean and standard deviation of intrinsic parameters by subset size.

input:
- results [subset-study results Dict]
- subset_sizes [images per subset]

outputs:
- Summary Dict [parameter name to subset sizes, means, and standard deviations]
"""
function summarise_by_subset(results::Dict{Int,Dict}, subset_sizes::Vector{Int})
    param_names = ["fx", "fy", "cx", "cy", "k1", "k2", "p1", "p2", "k3"]
    accessors   = [
        (K, d) -> K[1,1], (K, d) -> K[2,2],   # fx, fy
        (K, d) -> K[1,3], (K, d) -> K[2,3],   # cx, cy
        (K, d) -> d[1],   (K, d) -> d[2],     # k1, k2
        (K, d) -> d[3],   (K, d) -> d[4],     # p1, p2
        (K, d) -> d[5],                        # k3
    ]
    summary = Dict{String,NamedTuple}()
    for (name, acc) in zip(param_names, accessors)
        ns, μs, σs = Int[], Float64[], Float64[]
        for n in subset_sizes
            haskey(results, n) || continue
            d    = results[n]
            vals = [acc(K, dist)
                    for (K, dist) in zip(d["K_samples"], d["dist_samples"])
                    if length(dist) >= 5]
            isempty(vals) && continue
            push!(ns, n); push!(μs, mean(vals)); push!(σs, std(vals))
        end
        summary[name] = (n=ns, mean=μs, std=σs)
    end
    return summary
end

# =============================================================================
# ArUco corner uncertainty
# Subpixel edge fitting and per-corner covariance estimation.
# =============================================================================




# Functions for defining uncertainty in ArUco corner position 

# DEPRECATED: unused as of the uncertainty reorganization; candidate for deletion.
#=
function aruco_uncertainty_setup()
    # TODO: decide if I need this as its own function
    
    # OpenCV is 0 based; Julia is 1 based
    P = corners_px .+ 1
    # 1 = black, 0 = white; invert for intensity profile fitting
    Ig = 1 .- Float64.(Gray.(image_file));
    # interpolate intensity for sub-pixel sampling
    itp = linear_interpolation((1:size(Ig, 1), 1:size(Ig, 2)), 
							   Ig; 
							   extrapolation_bc = NaN
							  );  # itp(row, col), NaN outside
    ncols, nrows = 2, 2

    return P, Ig, itp, ncols, nrows
end
=#

function sigmoid(d::AbstractVector{<:Real}, p::AbstractVector{<:Real})
    """
    Sigmoid (logistic) edge model.

    Arguments:
    - d::AbstractVector{<:Real} - Positions (pixel offsets) at which the model is evaluated.
    - p::AbstractVector{<:Real} - Parameters [baseline, amplitude, edge location, width].
    
	(the amplitude is signed, so rising and falling edges are both covered)

    Returns:
    - AbstractVector - Model values at each position in d.
    """

    return p[1] .+ p[2] ./ (1 .+ exp.(-(d .- p[3]) ./ p[4]))

end

function fit_corner(intensity_interpolant, opencv_corner::AbstractVector{<:Real}, edge_direction_a::AbstractVector{<:Real}, edge_direction_b::AbstractVector{<:Real}; profile_positions::AbstractVector{<:Real}, half::Real)
    """
	Estimates where a corner is located using intensity profiles.
	-> Samples intensity profiles along several positions across a x and y edge 
	-> Fits a sigmoid to each profile; the sigmoid centre estimates where the edge is located
	-> Fits a straight line across those edge positions as a function of distance along the edge 
	-> Intersects the two fitted edge lines to estimate where the corner is 
	-> Compares the intensity-fitted corner coordinate to OpenCV's detected corner 
	-> Estimates a corner correction q relative to the OpenCV coordinates 
	-> Estimates covariance for the corner correction 

    Arguments:
    - intensity_interpolant - Image intensity interpolant, called as (row, col); should return NaN outside the image.
    - opencv_corner::AbstractVector{<:Real} - OpenCV corner [x, y] (Julia convention).
    - edge_direction_a::AbstractVector{<:Real} - Unit tangent vector [x] along edge A.
    - edge_direction_b::AbstractVector{<:Real} - Unit tangent vector [y] along edge B.
    - profile_positions::AbstractVector{<:Real} - Distances along each edge at which to sample profiles [px].
    - half::Real: Half-width of each profile across the edge [px].

    Returns:
	Nothing, if either edges has fewer than four usable profiles. Otherwise, a named tuple: 
	- q - Estimated edge-based corner minus opencv_corner
	- Σ - Estimated covariance matrix of q
	- edges - fitted information for each edge
	- d - offsets across sampled profiles 
 
    """

	# Offsets across each edge profile, centered on the approximate corner 
    profile_offsets = collect(-half:1.0:half)
    edge_fit_results = []

    for edge_tangent in (edge_direction_a, edge_direction_b)
		# a perpendicular unit normal to this edge 
        edge_normal = [-edge_tangent[2], edge_tangent[1]]
        usable_profiles = []

        for distance_along_edge in profile_positions
			# starts at this position along the edge; samples across normal
            profile_centre = opencv_corner .+ distance_along_edge .* edge_tangent
            profile_intensities = intensity_interpolant.(
				profile_centre[2] .+ profile_offsets .* edge_normal[2], 
				profile_centre[1] .+ profile_offsets .* edge_normal[1]
			)

			# skip profiles that would require sampling outside the image 
            any(isnan, profile_intensities) && continue

			# estimate the 2 plateaus and use the steepest step to initilise the sigmoid's baseline, signed amplitude, centre and width
            plateau_sample_count  = max(3, length(profile_intensities) ÷ 5)
            left_plateau = mean(profile_intensities[1:plateau_sample_count])
            right_plateau = mean(profile_intensities[end-plateau_sample_count+1:end])
            steepest_step_index  = argmax(abs.(diff(profile_intensities)))
			
            initial_parameters = [
				left_plateau, 
				right_plateau - left_plateau, 
				(profile_offsets[steepest_step_index] + profile_offsets[steepest_step_index + 1]) / 2, 
				1.0
			]

			# fit the sigmoid; constrain its centre + width to the sampled range 
            sigmoid_fit = try
                curve_fit(
					sigmoid, 
					profile_offsets, 
					profile_intensities, 
					initial_parameters; 
					lower = [-Inf, -Inf, -half, 0.3], 					
					upper = [Inf, Inf, half, 10.0]
				)
            catch
                nothing
            end

			# a failed fit or a v broad transition is not a useful edge
            (sigmoid_fit === nothing || sigmoid_fit.param[4] > 9.5) && continue 
			
            push!(usable_profiles, 
				  (distance_along_edge = distance_along_edge,
				   profile_intensities = profile_intensities, 
				   parameters = sigmoid_fit.param
				  )
				 )
        end

		# two fitted line parameters need more than two edge-position samples 
        length(usable_profiles) < 4 && return nothing

		# sigmoid centre (param 3) estimates edge offset for each profile
        fitted_edge_offsets = [
			profile.parameters[3] for profile in usable_profiles
		]
        fitted_distances = [profile.distance_along_edge for profile in usable_profiles
		]

		# model edge offset as intercept + slope * distance along edge 
        line_design_matrix = [ones(length(fitted_distances)) fitted_distances]
        # fit line coefficients (intercept, slope) to the edge offsets
        line_coefficients = line_design_matrix \ fitted_edge_offsets
        # compute residuals of the fitted line to estimate the scatter of the edge offsets
		line_residuals = fitted_edge_offsets .- line_design_matrix * line_coefficients

		# estimate line-coefficient covariance from residual scatter
        residual_variance = sum(abs2, line_residuals) / 
							(length(fitted_distances) - 2)
        line_coefficient_covariance = 
			residual_variance * inv(line_design_matrix' * line_design_matrix)
			
        push!(
			edge_fit_results, 
			  (
				  t = edge_tangent, 
				  n = edge_normal, 
				  a0 = line_coefficients[1],
				  b = line_coefficients[2], 
				  σ = sqrt(line_coefficient_covariance[1, 1]), 
                  residuals = line_residuals,
                  residual_variance = residual_variance,
				  profiles = usable_profiles
			  )
		)
    end

	first_edge = edge_fit_results[1]
	second_edge = edge_fit_results[2]

	# for a point q on an edge, its normal offset is a0 + b times its position along the tangent 
	# rearranging gives a linear equation per edge

	edge_intersection_matrix = [
		(first_edge.n .- first_edge.b .* first_edge.t)'; 
		(second_edge.n .- second_edge.b .* second_edge.t)'
	]
	edge_intercepts = [
		first_edge.a0, second_edge.a0
	]

	# solve the two edge equationsfor the correction from the OpenCV corner 
	corner_correction = edge_intersection_matrix \ edge_intercepts 

	# propagate the two intercept variances through the line intersection
	inverse_intersection_matrix = inv(edge_intersection_matrix)
	correction_covariance = inverse_intersection_matrix * Diagonal([first_edge.σ^2, second_edge.σ^2]) * inverse_intersection_matrix'

    return (q = corner_correction, 
			Σ = correction_covariance, 
			edges = edge_fit_results, 
			d = profile_offsets,
            opencv_corner = opencv_corner
		   )

end

function load_or_fit_aruco_corners(
    K, dist;
    extrinsic_dir, zoom_level, corner_results_file, aruco_image_file,
)
    cached_paths = (
        joinpath(extrinsic_dir, "aruco_corner_results_$(zoom_level).jld2"),
        corner_results_file,
        joinpath(
            extrinsic_dir,
            "extrinsic_bootstrap_$(zoom_level)_fit_corner_N1000.jld2",
        ),
    )
    for path in cached_paths
        isfile(path) || continue
        cached = JLD2.load(path)
        haskey(cached, "Σ_corners") || continue
        corners_key = haskey(cached, "fitted_corners_px") ? "fitted_corners_px" :
                      haskey(cached, "corners_px") ? "corners_px" : nothing
        corners_key === nothing && continue
        haskey(cached, "marker_ids") || path != corner_results_file ||
            error("$path is missing marker_ids.")
        @info "Loading cached ArUco corners and covariance" path
        return Matrix{Float64}(cached[corners_key]),
               Matrix{Float64}[Matrix{Float64}(Σ) for Σ in cached["Σ_corners"]]
    end

    isfile(aruco_image_file) || error("ArUco image not found: $aruco_image_file")
    corners_px, marker_ids = detect_aruco_corners(aruco_image_file, K, dist)
    size(corners_px) == (4, 2) ||
        error("Expected four ArUco corners, received $(size(corners_px)).")

    image = FileIO.load(aruco_image_file)
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
        result = fit_corner(
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
    mkpath(dirname(corner_results_file))
    JLD2.jldsave(
        corner_results_file;
        corner_results,
        opencv_corners_px=corners_px,
        fitted_corners_px=fitted_corners,
        Σ_corners,
        marker_ids,
        zoom_level,
    )
    @info "Saved regenerated ArUco observations and corner covariance" corner_results_file
    return fitted_corners, Σ_corners
end

# =============================================================================
# Extrinsic pose uncertainty
# Corner-noise pose bootstrap and pose statistics.
# =============================================================================

"""
Bootstrap ArUco marker extrinsics by perturbing detected corners and solving pose.
input:
- corners_px [4×2 matrix of detected image corners, px]
- K [3×3 camera matrix]
- dist [distortion coefficients]
- marker_length_mm [marker side length, mm]
- Σ_corners [vector of four 2×2 corner covariance matrices, px²]
- N [number of bootstrap samples]
- seed [random seed]
outputs:
- R_samples [rotation matrices, one 3×3 matrix per sample]
- t_samples [translation vectors, one 3-vector per sample, mm]
"""
function bootstrap_extrinsics(corners_px::Matrix{Float64}, K, dist,
                              marker_length_mm::Real;
                              Σ_corners::Vector{Matrix{Float64}},
                              N::Int, seed::Int)
    size(corners_px) == (4, 2) ||
        throw(ArgumentError("corners_px must be a 4×2 matrix."))
    marker_length_mm > 0 ||
        throw(ArgumentError("marker_length_mm must be positive."))
    length(Σ_corners) == 4 ||
        throw(ArgumentError("Σ_corners must contain exactly four covariance matrices."))
    for (i, Σ) in enumerate(Σ_corners)
        size(Σ) == (2, 2) ||
            throw(ArgumentError("Σ_corners[$i] must be a 2×2 matrix."))
        all(isfinite, Σ) ||
            throw(ArgumentError("Σ_corners[$i] must contain only finite values."))
        isapprox(Σ, transpose(Σ); atol=1e-12, rtol=1e-10) ||
            throw(ArgumentError("Σ_corners[$i] must be symmetric."))
        min_eigenvalue = (Σ[1, 1] + Σ[2, 2] -
                          sqrt((Σ[1, 1] - Σ[2, 2])^2 + 4Σ[1, 2]^2)) / 2
        min_eigenvalue >= -1e-10 * max(maximum(abs, Σ), 1.0) ||
            throw(ArgumentError("Σ_corners[$i] must be positive semidefinite."))
    end
    N >= 1 || throw(ArgumentError("N must be at least 1."))

    py"""
import cv2, numpy as np, time

corners = np.array($corners_px, dtype=np.float64).reshape(4, 2)
corner_covariances = np.asarray($Σ_corners, dtype=np.float64).reshape(4, 2, 2)
camera_matrix = np.array($K, dtype=np.float64).reshape(3, 3)
dist_coeffs = np.array($dist, dtype=np.float64).reshape(-1)
marker_length = float($marker_length_mm)
n_samples = int($N)
rng = np.random.default_rng(int($seed))
obj_points = np.array([
    [-marker_length/2,  marker_length/2, 0],
    [ marker_length/2,  marker_length/2, 0],
    [ marker_length/2, -marker_length/2, 0],
    [-marker_length/2, -marker_length/2, 0],
], dtype=np.float64)

R_samples, t_samples = [], []
print(f"Starting bootstrap: N={n_samples} samples...", flush=True)
start_time = time.perf_counter()
progress_interval = max(1, min(100, round(n_samples / 10)))
for i in range(n_samples):
    noisy_corners = corners.copy()
    for corner_index in range(4):
        noisy_corners[corner_index] += rng.multivariate_normal(
            np.zeros(2), corner_covariances[corner_index], check_valid="raise")
    success, rvec, tvec = cv2.solvePnP(
        obj_points, noisy_corners, camera_matrix, dist_coeffs)
    if not success:
        raise RuntimeError(f"solvePnP failed for bootstrap sample {i + 1}")
    R, _ = cv2.Rodrigues(rvec)
    R_samples.append(R.tolist())
    t_samples.append(tvec.reshape(-1).tolist())
    samples_done = i + 1
    if samples_done % progress_interval == 0 or samples_done == n_samples:
        elapsed = time.perf_counter() - start_time
        samples_remaining = n_samples - samples_done
        eta = elapsed / samples_done * samples_remaining
        percent = round(100 * samples_done / n_samples)
        print(
            f"  sample {samples_done}/{n_samples}  ({percent}%)  "
            f"elapsed: {elapsed:.1f} s  ETA: {eta:.1f} s",
            flush=True)
total_elapsed = time.perf_counter() - start_time
print(f"Bootstrap complete: {n_samples} samples in {total_elapsed:.1f} s", flush=True)
result = (R_samples, t_samples)
"""
    R_raw, t_raw = py"result"
    R_samples = Matrix{Float64}[
        Matrix{Float64}(R_raw[i, :, :]) for i in axes(R_raw, 1)
    ]
    t_samples = Vector{Float64}[
        Float64.(vec(t_raw[i, :])) for i in axes(t_raw, 1)
    ]
    return R_samples, t_samples
end


"""
Compute sample means and covariance matrices for extrinsic pose samples.
input:
- R_samples [rotation matrices, one 3×3 matrix per sample]
- t_samples [translation vectors, one 3-vector per sample, mm]
outputs:
- Σ_t [3×3 translation covariance matrix, mm²]
- Σ_rot [3×3 covariance matrix of Δrvec relative to R_mean, rad²]
- rvec_mean [mean rotation vector, rad]
- t_mean [mean translation vector, mm]
- θ_samples [geodesic rotation distances from the mean, degrees]
- σ_rotation [standard deviation of geodesic rotation distance, degrees]
- R_mean [mean rotation matrix, 3×3]
"""
function extrinsic_covariance(R_samples::Vector{Matrix{Float64}},
                              t_samples::Vector{Vector{Float64}})
    length(R_samples) == length(t_samples) ||
        throw(ArgumentError("R_samples and t_samples must have the same number of samples."))
    length(R_samples) >= 2 ||
        throw(ArgumentError("At least two pose samples are required to estimate covariance."))
    all(size(R) == (3, 3) for R in R_samples) ||
        throw(ArgumentError("Every rotation sample must be a 3×3 matrix."))
    all(length(t) == 3 for t in t_samples) ||
        throw(ArgumentError("Every translation sample must contain three values."))

    py"""
import cv2, numpy as np

rotations = $R_samples
rvec_samples = []
for rotation in rotations:
    rotation_matrix = np.array(rotation, dtype=np.float64)
    rvec, _ = cv2.Rodrigues(rotation_matrix)
    rvec_samples.append(rvec.reshape(-1).tolist())
"""
    rvec_raw = py"rvec_samples"
    rvec_matrix = Matrix{Float64}(rvec_raw)
    t_matrix = reduce(vcat, permutedims.(t_samples))
    rvec_mean = vec(mean(rvec_matrix; dims=1))
    t_mean = vec(mean(t_matrix; dims=1))
    Σ_t = Matrix{Float64}(cov(t_matrix; dims=1))

    py"""
import cv2, numpy as np

rotations = $R_samples
mean_rvec = np.array($rvec_mean, dtype=np.float64).reshape(3, 1)
mean_rotation, _ = cv2.Rodrigues(mean_rvec)
delta_rvec_samples = []
for rotation in rotations:
    delta_rotation = mean_rotation.T @ np.array(rotation, dtype=np.float64)
    delta_rvec, _ = cv2.Rodrigues(delta_rotation)
    delta_rvec_samples.append(delta_rvec.reshape(-1).tolist())
result = (mean_rotation.tolist(), delta_rvec_samples)
"""
    R_mean_raw, delta_rvec_raw = py"result"
    R_mean = Matrix{Float64}(R_mean_raw)
    delta_rvec_matrix = Matrix{Float64}(delta_rvec_raw)
    Σ_rot = Matrix{Float64}(cov(delta_rvec_matrix; dims=1))
    θ_samples = Float64[]
    sizehint!(θ_samples, length(R_samples))
    for R in R_samples
        ΔR = transpose(R_mean) * R
        cosine = clamp((ΔR[1, 1] + ΔR[2, 2] + ΔR[3, 3] - 1) / 2, -1.0, 1.0)
        push!(θ_samples, acosd(cosine))
    end
    σ_rotation = std(θ_samples)
    return Σ_t, Σ_rot, rvec_mean, t_mean, θ_samples, σ_rotation, R_mean
end
"""
Apply a Rodrigues rotation-vector increment as a right perturbation of `R_mean`.

The increment is expressed in the mean camera rotation's local frame.
"""
function _right_perturbed_rotation(R_mean::AbstractMatrix, delta::AbstractVector)
    length(delta) == 3 || throw(ArgumentError("Rotation increment must have length 3."))
    θ2 = sum(abs2, delta)
    T = typeof(θ2)
    if θ2 < 1e-12
        a = one(θ2) - θ2 / 6 + θ2^2 / 120
        b = one(θ2) / 2 - θ2 / 24 + θ2^2 / 720
    else
        θ = sqrt(θ2)
        a = sin(θ) / θ
        b = (one(θ) - cos(θ)) / θ2
    end
    W = T[
        zero(θ2) -delta[3] delta[2]
        delta[3] zero(θ2) -delta[1]
        -delta[2] delta[1] zero(θ2)
    ]
    I3 = Matrix{T}(I, 3, 3)
    return R_mean * (I3 + a .* W + b .* (W * W))
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
        R_mean = _right_perturbed_rotation(R_mean, mean_increment)
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

# =============================================================================
# Camera-to-drop projection
# Forward and inverse mapping between image pixels and the model plane.
# =============================================================================


"""Extract a real scalar value to anchor an implicit undistortion solve."""
_camtodrop_primal(x::ForwardDiff.Dual) = _camtodrop_primal(ForwardDiff.value(x))
_camtodrop_primal(x::Real) = Float64(x)


"""
Return the model-local (r,z) coordinates for fixed pixels.

`parameters` is `[fx, fy, cx, cy, k1, k2, p1, p2, k3,
δrx, δry, δrz, δtx, δty, δtz, δoffset]`. Rotation increments are right
perturbations of `R_mean`, matching `extrinsic_covariance`'s
`R_mean' * R_sample` convention. The relative-pose offset shifts `t_rel` along
the model plane normal. For dual-number inputs, distortion inversion uses the
implicit-function derivative evaluated at `reference_parameters`; it does not
differentiate through the iterative undistortion solver. `reference_pixels`
anchors the normalized-ray solution when evaluating perturbed pixels.
"""
function inverse_project_camtodrop(
    u::AbstractVector,
    v::AbstractVector,
    parameters::AbstractVector,
    R_mean::AbstractMatrix,
    t_mean::AbstractVector,
    R_rel::AbstractMatrix,
    t_rel::AbstractVector;
    plane_coordinate::Real = 0.0,
    reference_parameters::Union{Nothing,AbstractVector} = nothing,
    reference_pixels::Union{Nothing,AbstractVector} = nothing,
    iterations::Int = 20,
)
    length(u) == length(v) ||
        throw(ArgumentError("u and v must have the same length."))
    length(parameters) == 16 ||
        throw(ArgumentError("parameters must contain 16 values."))
    length(t_mean) == 3 && length(t_rel) == 3 ||
        throw(ArgumentError("Camera and relative translations must have length 3."))
    iterations >= 1 || throw(ArgumentError("iterations must be positive."))

    is_float_input =
        eltype(parameters) <: AbstractFloat &&
        eltype(u) <: AbstractFloat &&
        eltype(v) <: AbstractFloat
    if is_float_input
        fx, fy, cx, cy = parameters[1:4]
        K = [fx 0.0 cx; 0.0 fy cy; 0.0 0.0 1.0]
        dist = Float64.(parameters[5:9])
        R = _right_perturbed_rotation(
            Float64.(R_mean), Float64.(parameters[10:12]),
        )
        t = Float64.(t_mean) .+ Float64.(parameters[13:15])
        R_rel_float = Float64.(R_rel)
        t_rel_offset = Float64.(t_rel) .+
                       Float64(parameters[16]) .* R_rel_float[:, 3]
        r, z, _ = Camera.pixel_to_model(
            u, v, R, t, K, dist, R_rel_float, t_rel_offset;
            plane_coordinate=plane_coordinate,
            iterations=iterations,
        )
        return r, z
    end

    reference_parameters === nothing && throw(ArgumentError(
        "reference_parameters is required for ForwardDiff inputs.",
    ))
    length(reference_parameters) == 16 ||
        throw(ArgumentError("reference_parameters must contain 16 values."))
    reference_pixels === nothing || length(reference_pixels) == 2length(u) ||
        throw(ArgumentError("reference_pixels must contain interleaved u,v values."))

    p0 = Float64.(reference_parameters)
    fx0, fy0, cx0, cy0 = p0[1:4]
    dist0 = p0[5:9]
    R = _right_perturbed_rotation(R_mean, view(parameters, 10:12))
    t = t_mean .+ view(parameters, 13:15)
    T = promote_type(eltype(parameters), eltype(u), eltype(v))
    r = Vector{T}(undef, length(u))
    z = similar(r)
    normal = Float64.(R_rel[:, 3])
    t_relative = t_rel .+ parameters[16] .* normal
    plane_point = t_relative .+ plane_coordinate .* normal
    camera_center = -(transpose(R) * t)

    for i in eachindex(u, v)
        reference_u = reference_pixels === nothing ?
            _camtodrop_primal(u[i]) : reference_pixels[2i - 1]
        reference_v = reference_pixels === nothing ?
            _camtodrop_primal(v[i]) : reference_pixels[2i]
        xd0 = (reference_u - cx0) / fx0
        yd0 = (reference_v - cy0) / fy0
        x0_vec, y0_vec = Camera.undistortion(
            [xd0], [yd0], dist0; max_iter=iterations, tol=1e-12,
        )
        x0, y0 = x0_vec[1], y0_vec[1]
        Jx = ForwardDiff.jacobian(
            q -> collect(Camera.distortion(q[1], q[2], dist0)),
            Float64[x0, y0],
        )
        abs(det(Jx)) > eps(Float64) ||
            throw(ArgumentError("Distortion Jacobian is singular at pixel $i."))

        xd = (u[i] - parameters[3]) / parameters[1]
        yd = (v[i] - parameters[4]) / parameters[2]
        d0 = Camera.distortion(x0, y0, dist0)
        d = Camera.distortion(x0, y0, view(parameters, 5:9))
        residual = [
            d[1] - xd - (d0[1] - xd0),
            d[2] - yd - (d0[2] - yd0),
        ]
        x = x0 - (Jx \ residual)[1]
        y = y0 - (Jx \ residual)[2]

        ray_camera = [x, y, one(x)]
        ray_marker = transpose(R) * ray_camera
        denominator = dot(normal, ray_marker)
        abs(denominator) > 1e-12 ||
            throw(ArgumentError("Viewing ray is parallel to the droplet plane at pixel $i."))
        λ = dot(normal, plane_point .- camera_center) / denominator
        point_marker = camera_center .+ λ .* ray_marker
        point_model = transpose(R_rel) * (point_marker .- t_relative)
        r[i] = point_model[1]
        z[i] = point_model[2]
    end
    return r, z
end

"""
Forward-project model-local contour coordinates using the same transform order
as `Model.ForwardModel`. The meridian coordinate is fixed by `plane_coordinate`.
Returned pixels use the OpenCV zero-based coordinate convention of the intrinsics.
"""
function forward_project_camtodrop(
    r::AbstractVector,
    z::AbstractVector,
    parameters::AbstractVector,
    R_mean::AbstractMatrix,
    t_mean::AbstractVector,
    R_rel::AbstractMatrix,
    t_rel::AbstractVector;
    plane_coordinate::Real = 0.0,
)
    length(r) == length(z) ||
        throw(ArgumentError("r and z must have the same length."))
    length(parameters) == 16 ||
        throw(ArgumentError("parameters must contain 16 values."))

    fx, fy, cx, cy = parameters[1:4]
    K = [fx 0.0 cx; 0.0 fy cy; 0.0 0.0 1.0]
    dist = parameters[5:9]
    R = _right_perturbed_rotation(R_mean, view(parameters, 10:12))
    t = t_mean .+ view(parameters, 13:15)
    t_relative = t_rel .+ parameters[16] .* R_rel[:, 3]

    if eltype(parameters) <: AbstractFloat &&
       eltype(r) <: AbstractFloat && eltype(z) <: AbstractFloat
        return Camera.project_model_to_image(
            Float64.(r), Float64.(z),
            fill(Float64(plane_coordinate), length(r)),
            Float64.(R), Float64.(t), Float64.(K), Float64.(dist),
            Float64.(R_rel), Float64.(t_relative),
        )
    end

    R_total = R * R_rel
    t_total = R * t_relative + t
    T = promote_type(eltype(r), eltype(z), eltype(parameters))
    X = vcat(
        reshape(T.(r), 1, :),
        reshape(T.(z), 1, :),
        fill(one(parameters[1]) * plane_coordinate, 1, length(r)),
    )
    X_camera = R_total * X .+ reshape(t_total, 3, 1)
    x = X_camera[1, :] ./ X_camera[3, :]
    y = X_camera[2, :] ./ X_camera[3, :]
    xd, yd = Camera.distortion(x, y, dist)
    return fx .* xd .+ cx, fy .* yd .+ cy, X_camera[3, :]
end

# =============================================================================
# Uncertainty propagation and sampling
# Projection Jacobians, covariance assembly, and uncertainty sampling.
# =============================================================================

"""
Calculate inverse-projection Jacobians with respect to the nine intrinsics,
six extrinsic parameters, the relative-pose normal offset, and the pixel data.
Undistortion derivatives use the implicit function theorem.
"""
function camtodrop_jacobians(
    u::AbstractVector,
    v::AbstractVector,
    intrinsics::AbstractVector,
    R_mean::AbstractMatrix,
    t_mean::AbstractVector,
    R_rel::AbstractMatrix,
    t_rel::AbstractVector;
    plane_coordinate::Real = 0.0,
)
    length(intrinsics) == 9 ||
        throw(ArgumentError("intrinsics must follow [fx,fy,cx,cy,k1,k2,p1,p2,k3]."))
    base = vcat(Float64.(intrinsics), zeros(6), 0.0)
    reference = copy(base)
    n = length(u)

    flatten_rz(r, z) = reduce(vcat, ([r[i], z[i]] for i in eachindex(r, z)))
    inverse_parameters(p, pixels) = begin
        r, z = inverse_project_camtodrop(
            view(pixels, 1:2:2n), view(pixels, 2:2:2n), p,
            R_mean, t_mean, R_rel, t_rel;
            plane_coordinate=plane_coordinate,
            reference_parameters=reference,
            reference_pixels=pixels0,
        )
        flatten_rz(r, z)
    end
    pixels0 = reduce(vcat, ([u[i], v[i]] for i in eachindex(u, v)))
    y0 = inverse_parameters(base, pixels0)

    J_intrinsics = ForwardDiff.jacobian(
        intrinsic_sample ->
            inverse_parameters(vcat(intrinsic_sample, zeros(6), 0.0), pixels0),
        Float64.(intrinsics),
    )
    J_extrinsics = ForwardDiff.jacobian(
        extrinsic_sample ->
            inverse_parameters(vcat(Float64.(intrinsics), extrinsic_sample, 0.0), pixels0),
        zeros(6),
    )
    J_offset = ForwardDiff.jacobian(
        offset_sample ->
            inverse_parameters(vcat(Float64.(intrinsics), zeros(6), offset_sample), pixels0),
        [0.0],
    )
    J_pixels = ForwardDiff.jacobian(
        pixel_sample -> inverse_parameters(base, pixel_sample),
        pixels0,
    )

    return (coordinates=y0, J_intrinsics=J_intrinsics,
            J_extrinsics=J_extrinsics, J_offset=J_offset,
            J_pixels=J_pixels, parameter_order=(
                "fx", "fy", "cx", "cy", "k1", "k2", "p1", "p2", "k3",
                "δrx_right_rad", "δry_right_rad", "δrz_right_rad",
                "δtx_mm", "δty_mm", "δtz_mm", "δoffset_normal_mm",
            ))
end

"""
Return Jacobian of undistorted normalized coordinates with respect to the
five OpenCV distortion coefficients, using `∂x/∂φ = -(∂d/∂x)⁻¹∂d/∂φ`.
"""
function undistortion_implicit_parameter_jacobian(
    xd::Real,
    yd::Real,
    dist_coeffs::AbstractVector{<:Real},
)
    length(dist_coeffs) >= 5 ||
        throw(ArgumentError("At least five OpenCV distortion coefficients are required."))
    dist = Float64.(dist_coeffs[1:5])
    x_values, y_values = Camera.undistortion(
        [Float64(xd)], [Float64(yd)], dist; max_iter=40, tol=1e-13,
    )
    x0 = Float64[x_values[1], y_values[1]]
    Jx = ForwardDiff.jacobian(
        x -> collect(Camera.distortion(x[1], x[2], dist)), x0,
    )
    Jφ = ForwardDiff.jacobian(
        coefficients -> collect(Camera.distortion(x0[1], x0[2], coefficients)),
        dist,
    )
    return -(Jx \ Jφ)
end

"""
Sample a multivariate normal using a positive-semidefinite covariance matrix.
"""
function sample_gaussian(rng::AbstractRNG, mean_vector::AbstractVector,
                         covariance::AbstractMatrix)
    n = length(mean_vector)
    size(covariance) == (n, n) ||
        throw(ArgumentError("Covariance shape does not match the mean vector."))
    isapprox(covariance, transpose(covariance); atol=1e-10, rtol=1e-8) ||
        throw(ArgumentError("Covariance matrix must be symmetric."))
    decomposition = eigen(Symmetric(Float64.(covariance)))
    tolerance = 1e-10 * max(maximum(abs, decomposition.values), 1.0)
    minimum(decomposition.values) >= -tolerance ||
        throw(ArgumentError("Covariance matrix must be positive semidefinite."))
    factor = decomposition.vectors * Diagonal(sqrt.(max.(decomposition.values, 0.0)))
    return Float64.(mean_vector) .+ factor * randn(rng, n)
end

"""
Run one solvePnP bootstrap draw per intrinsic-parameter sample, matching
`bootstrap_extrinsics`'s target geometry, corner perturbations, and solver.
"""
function joint_intrinsic_extrinsic_bootstrap(
    corners_px::AbstractMatrix{<:Real},
    corner_covariances::AbstractVector{<:AbstractMatrix},
    intrinsic_samples::AbstractMatrix{<:Real},
    marker_length_mm::Real;
    seed::Int,
)
    size(corners_px) == (4, 2) ||
        throw(ArgumentError("corners_px must be a 4×2 matrix."))
    length(corner_covariances) == 4 ||
        throw(ArgumentError("Exactly four corner covariance matrices are required."))
    size(intrinsic_samples, 2) == 9 ||
        throw(ArgumentError("intrinsic_samples must have nine columns."))
    marker_length_mm > 0 ||
        throw(ArgumentError("marker_length_mm must be positive."))

    py"""
import cv2, numpy as np

corners = np.asarray($corners_px, dtype=np.float64).reshape(4, 2)
corner_covariances = np.asarray($corner_covariances, dtype=np.float64).reshape(4, 2, 2)
samples = np.asarray($intrinsic_samples, dtype=np.float64)
marker_length = float($marker_length_mm)
rng = np.random.default_rng(int($seed))
obj_points = np.array([
    [-marker_length/2,  marker_length/2, 0],
    [ marker_length/2,  marker_length/2, 0],
    [ marker_length/2, -marker_length/2, 0],
    [-marker_length/2, -marker_length/2, 0],
], dtype=np.float64)
rotations, translations = [], []
for i, phi in enumerate(samples):
    K = np.array([
        [phi[0], 0.0, phi[2]],
        [0.0, phi[1], phi[3]],
        [0.0, 0.0, 1.0],
    ], dtype=np.float64)
    dist = np.asarray(phi[4:9], dtype=np.float64)
    noisy_corners = corners.copy()
    for corner_index in range(4):
        noisy_corners[corner_index] += rng.multivariate_normal(
            np.zeros(2), corner_covariances[corner_index], check_valid="raise")
    success, rvec, tvec = cv2.solvePnP(obj_points, noisy_corners, K, dist)
    if not success:
        raise RuntimeError(f"solvePnP failed for joint bootstrap sample {i + 1}")
    R, _ = cv2.Rodrigues(rvec)
    rotations.append(R.tolist())
    translations.append(tvec.reshape(-1).tolist())
result = (rotations, translations)
"""
    R_raw, t_raw = py"result"
    R_samples = Matrix{Float64}[
        Matrix{Float64}(R_raw[i, :, :]) for i in axes(R_raw, 1)
    ]
    t_samples = Vector{Float64}[
        Float64.(vec(t_raw[i, :])) for i in axes(t_raw, 1)
    ]
    return R_samples, t_samples
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

"""
Build the interleaved pixel covariance for correlated normal-only edge errors.
"""
function edge_normal_pixel_covariance(
    u::AbstractVector{<:Real},
    v::AbstractVector{<:Real},
    sigma_n_px::AbstractVector{<:Real},
    correlation::AbstractMatrix{<:Real},
)
    n = length(u)
    length(v) == n && length(sigma_n_px) == n ||
        throw(ArgumentError("u, v, and sigma_n_px lengths must match."))
    n >= 2 || throw(ArgumentError("At least two edge points are required."))
    all(isfinite, u) && all(isfinite, v) &&
        all(x -> isfinite(x) && x > 0, sigma_n_px) ||
        throw(ArgumentError("Edge coordinates must be finite and σ must be positive."))
    corr = _validate_correlation_matrix(correlation, n)

    normals = Matrix{Float64}(undef, n, 2)
    for i in 1:n
        left, right = i == 1 ? (1, 2) : i == n ? (n - 1, n) : (i - 1, i + 1)
        du = Float64(u[right] - u[left])
        dv = Float64(v[right] - v[left])
        magnitude = hypot(du, dv)
        magnitude > 0 ||
            throw(ArgumentError("Degenerate image-space tangent at edge point $i."))
        normals[i, :] .= (-dv / magnitude, du / magnitude)
    end

    covariance = zeros(Float64, 2n, 2n)
    for i in 1:n, j in 1:n
        block = (
            Float64(sigma_n_px[i]) * Float64(sigma_n_px[j]) * corr[i, j]
        ) .* (view(normals, i, :) * transpose(view(normals, j, :)))
        covariance[2i-1:2i, 2j-1:2j] .= block
    end
    return Matrix(Symmetric(covariance))
end

normal_sigmas(covariance, normals) =
    sqrt.(max.(diag(normal_covariance(covariance, normals)), 0.0))


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
    R_rel, t_rel; plane_coordinate::Real,
)
    parameters = vcat(phi, extrinsic_delta, pose_offset)
    r, z = inverse_project_camtodrop(
        pixels_u, pixels_v, parameters, R_mean, t_mean, R_rel, t_rel;
        plane_coordinate,
    )
    return stack_rz(r, z)
end

# =============================================================================
# Shape and position analysis
# Contour-shape metrics, shared-corner sensitivity, and relative-position effects.
# =============================================================================

"""
Check that inverse-projected points run along one meridian from apex to needle.
"""
function validate_observed_contour_order(
    r::AbstractVector{<:Real},
    z::AbstractVector{<:Real},
    needle_radius_mm::Real,
)
    length(r) == length(z) ||
        throw(ArgumentError("r and z lengths differ."))
    length(r) >= 4 ||
        throw(ArgumentError("At least four observed contour points are required."))
    needle_radius_mm > 0 && isfinite(needle_radius_mm) ||
        throw(ArgumentError("needle_radius_mm must be finite and positive."))
    all(isfinite, r) && all(isfinite, z) ||
        throw(ArgumentError("Observed model-plane contour must be finite."))

    tolerance = 0.05 * needle_radius_mm
    one_meridian_side =
        minimum(r) >= -tolerance || maximum(r) <= tolerance
    one_meridian_side ||
        throw(ArgumentError("Observed points must follow one meridian side."))
    abs(r[1]) <= 0.15 * needle_radius_mm ||
        throw(ArgumentError("Contour ordering must start near the apex."))
    abs(abs(r[end]) - needle_radius_mm) <= 0.15 * needle_radius_mm ||
        throw(ArgumentError("Contour ordering must end near the needle."))
    z[end] > z[1] ||
        throw(ArgumentError("Contour ordering must proceed from apex to needle."))
    return nothing
end

function apex_radius_measure(r::AbstractVector{<:Real}, z::AbstractVector{<:Real})
    length(r) == length(z) ||
        throw(ArgumentError("r and z lengths differ."))
    indices = sortperm(abs.(r))[1:min(16, length(r))]
    length(indices) >= 4 || return NaN
    design = hcat(ones(length(indices)), Float64.(r[indices]).^2)
    coefficient = design \ Float64.(z[indices])
    curvature = abs(coefficient[2])
    return curvature > eps(Float64) ? 1 / (2curvature) : Inf
end

"""
Measure the deterministic pose and contour effects of shared 0.1 px corner shifts.

The three cases apply the same shift to all corners in x, in y, or radially
outward from their centroid. The contour displacements are median absolute
normal components relative to the pose solved from the unshifted corners.
"""
function shared_corner_shift_sensitivity(
    corners_px::AbstractMatrix{<:Real},
    K::AbstractMatrix{<:Real},
    dist::AbstractVector{<:Real},
    marker_length_mm::Real,
    virtual_u::AbstractVector{<:Real},
    virtual_v::AbstractVector{<:Real},
    intrinsics::AbstractVector{<:Real},
    R_rel::AbstractMatrix{<:Real},
    t_rel::AbstractVector{<:Real},
    joint_camera_normal_sigma::AbstractVector{<:Real};
    plane_coordinate::Real = 0.0,
)
    size(corners_px) == (4, 2) ||
        throw(ArgumentError("corners_px must be a 4×2 matrix."))
    size(K) == (3, 3) ||
        throw(ArgumentError("K must be a 3×3 matrix."))
    length(intrinsics) == 9 ||
        throw(ArgumentError("intrinsics must contain nine values."))
    length(virtual_u) == length(virtual_v) ||
        throw(ArgumentError("virtual_u and virtual_v lengths differ."))
    length(virtual_u) >= 4 ||
        throw(ArgumentError("At least four virtual contour points are required."))
    size(R_rel) == (3, 3) && length(t_rel) == 3 ||
        throw(ArgumentError("Relative pose must contain a 3×3 rotation and 3-vector translation."))
    length(joint_camera_normal_sigma) == length(virtual_u) ||
        throw(ArgumentError("Joint camera sigma length must match the contour."))
    marker_length_mm > 0 ||
        throw(ArgumentError("marker_length_mm must be positive."))
    all(isfinite, corners_px) && all(isfinite, K) && all(isfinite, dist) &&
        all(isfinite, intrinsics) && all(isfinite, virtual_u) &&
        all(isfinite, virtual_v) && all(isfinite, joint_camera_normal_sigma) ||
        throw(ArgumentError("Inputs must contain only finite values."))
    median_joint_camera_sigma_mm = median(joint_camera_normal_sigma)
    median_joint_camera_sigma_mm > 0 ||
        throw(ArgumentError("Median joint camera sigma must be positive."))

    py"""
import cv2, numpy as np

corners = np.asarray($corners_px, dtype=np.float64).reshape(4, 2)
camera_matrix = np.asarray($K, dtype=np.float64).reshape(3, 3)
dist_coeffs = np.asarray($dist, dtype=np.float64).reshape(-1)
marker_length = float($marker_length_mm)
radial = corners - np.mean(corners, axis=0)
radial_norms = np.linalg.norm(radial, axis=1)
if np.any(radial_norms <= np.finfo(np.float64).eps):
    raise ValueError("A corner coincides with the marker centroid.")
outward = 0.1 * radial / radial_norms[:, None]
corner_shifts = [
    np.column_stack((np.full(4, 0.1), np.zeros(4))),
    np.column_stack((np.zeros(4), np.full(4, 0.1))),
    outward,
]
obj_points = np.array([
    [-marker_length/2,  marker_length/2, 0],
    [ marker_length/2,  marker_length/2, 0],
    [ marker_length/2, -marker_length/2, 0],
    [-marker_length/2, -marker_length/2, 0],
], dtype=np.float64)
rotations, translations = [], []
for shift in [np.zeros((4, 2))] + corner_shifts:
    success, rvec, tvec = cv2.solvePnP(
        obj_points, corners + shift, camera_matrix, dist_coeffs)
    if not success:
        raise RuntimeError("solvePnP failed for a shared-corner shift case.")
    rotation, _ = cv2.Rodrigues(rvec)
    rotations.append(rotation.tolist())
    translations.append(tvec.reshape(-1).tolist())
result = (rotations, translations)
"""
    R_raw, t_raw = py"result"
    rotations = Matrix{Float64}[
        Matrix{Float64}(R_raw[i, :, :]) for i in axes(R_raw, 1)
    ]
    translations = Vector{Float64}[
        Float64.(vec(t_raw[i, :])) for i in axes(t_raw, 1)
    ]
    R_reference, t_reference = rotations[1], translations[1]
    reference_parameters = vcat(Float64.(intrinsics), zeros(6), 0.0)
    r_reference, z_reference = inverse_project_camtodrop(
        virtual_u, virtual_v, reference_parameters,
        R_reference, t_reference, R_rel, t_rel;
        plane_coordinate=plane_coordinate,
    )
    normals = contour_normals(r_reference, z_reference)
    reference_apex_radius = apex_radius_measure(r_reference, z_reference)
    labels = (
        "a: all corners +0.1 px in x",
        "b: all corners +0.1 px in y",
        "c: all corners 0.1 px outward",
    )
    cases = map(2:4) do pose_index
        rotation_change = 1_000 .* rotation_log(
            transpose(R_reference) * rotations[pose_index],
        )
        translation_change = 1_000 .* (
            translations[pose_index] .- t_reference
        )
        contour = sample_inverse_contour(
            virtual_u, virtual_v, intrinsics,
            vcat(rotation_change ./ 1_000, translation_change ./ 1_000),
            0.0, R_reference, t_reference, R_rel, t_rel;
            plane_coordinate=plane_coordinate,
        )
        r_case, z_case = unstack_rz(contour)
        normal_displacement = (r_case .- r_reference) .* normals[:, 1] .+
                              (z_case .- z_reference) .* normals[:, 2]
        apex_radius_change =
            apex_radius_measure(r_case, z_case) - reference_apex_radius
        (
            label=labels[pose_index - 1],
            rotation_change_mrad=rotation_change,
            translation_change_um=translation_change,
            median_abs_normal_displacement_um=1_000 * median(abs.(normal_displacement)),
            apex_radius_change_um=1_000 * apex_radius_change,
        )
    end
    return (
        cases=cases,
        median_joint_camera_sigma_um=1_000 * median_joint_camera_sigma_mm,
    )
end
function relative_position_assessment(
    J_offset::AbstractMatrix,
    normals::AbstractMatrix,
    joint_camera_normal_sigma::AbstractVector,
    t_mean::AbstractVector,
)
    n = size(normals, 1)
    n > 0 || throw(ArgumentError("At least one contour normal is required."))
    size(normals, 2) == 2 ||
        throw(ArgumentError("Contour normals must have two columns."))
    size(J_offset) == (2n, 1) ||
        throw(ArgumentError("Offset Jacobian must have one column per contour."))
    length(joint_camera_normal_sigma) == n ||
        throw(ArgumentError("Joint-camera normal σ does not match the contour length."))
    length(t_mean) == 3 ||
        throw(ArgumentError("Mean translation must have three coordinates."))
    all(isfinite, J_offset) ||
        throw(ArgumentError("Offset Jacobian must contain finite values."))
    all(isfinite, normals) ||
        throw(ArgumentError("Contour normals must contain finite values."))
    all(x -> isfinite(x) && x >= 0, joint_camera_normal_sigma) ||
        throw(ArgumentError("Joint-camera normal σ must be finite and nonnegative."))
    all(isfinite, t_mean) ||
        throw(ArgumentError("Mean translation must contain finite values."))

    camera_marker_distance_mm = abs(Float64(t_mean[3]))
    camera_marker_distance_mm > 0 ||
        throw(ArgumentError("Camera-marker optical depth must be positive."))

    unit_offset_covariance = J_offset * transpose(J_offset)
    median_normal_sigma_per_mm_offset =
        median(normal_sigmas(unit_offset_covariance, normals))
    median_normal_sigma_per_mm_offset > 0 ||
        throw(ArgumentError("Offset Jacobian has zero median normal sensitivity."))
    median_joint_camera_normal_sigma_mm =
        median(Float64.(joint_camera_normal_sigma))
    break_even_offset_mm =
        median_joint_camera_normal_sigma_mm / median_normal_sigma_per_mm_offset

    return (
        camera_marker_distance_mm=camera_marker_distance_mm,
        median_normal_sigma_per_mm_offset=median_normal_sigma_per_mm_offset,
        median_joint_camera_normal_sigma_mm=median_joint_camera_normal_sigma_mm,
        break_even_offset_mm=break_even_offset_mm,
        relative_scale_error_per_mm=inv(camera_marker_distance_mm),
        relative_surface_tension_error_per_mm=2 / camera_marker_distance_mm,
    )
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

# =============================================================================
# Not uncertainty-related
# Detection, calibration, I/O, and plotting support outside uncertainty calculations.
# =============================================================================

"""
Load and validate 1-based edge pixels from an image with the same dimensions
as the raw image used for the ArUco pose. Relative `image_file` paths are
resolved from `image_root`.
"""
function load_edge_uncertainty_data(
    path::AbstractString;
    image_root::AbstractString,
    expected_image_file::AbstractString,
)
    isfile(path) || throw(ArgumentError("Edge-data file not found: $path"))
    data = JLD2.load(path)
    require_keys(
        data, ("u", "v", "sigma_n_px", "sigma_description", "image_file"), path,
    )

    u, v, sigma_n_px = data["u"], data["v"], data["sigma_n_px"]
    u isa Vector{Float64} ||
        throw(ArgumentError("Edge-data key \"u\" must be Vector{Float64}."))
    v isa Vector{Float64} ||
        throw(ArgumentError("Edge-data key \"v\" must be Vector{Float64}."))
    sigma_n_px isa Vector{Float64} ||
        throw(ArgumentError("Edge-data key \"sigma_n_px\" must be Vector{Float64}."))
    sigma_description = data["sigma_description"]
    sigma_description isa String ||
        throw(ArgumentError("Edge-data key \"sigma_description\" must be a String."))
    !isempty(strip(sigma_description)) ||
        throw(ArgumentError("sigma_description must be non-empty."))
    image_file = data["image_file"]
    image_file isa String ||
        throw(ArgumentError("Edge-data key \"image_file\" must be a String."))

    n = length(u)
    n >= 4 || throw(ArgumentError("At least four edge points are required."))
    length(v) == n && length(sigma_n_px) == n ||
        throw(ArgumentError("u, v, and sigma_n_px lengths must match."))
    all(isfinite, u) && all(isfinite, v) ||
        throw(ArgumentError("Edge coordinates must be finite."))
    all(x -> isfinite(x) && x > 0, sigma_n_px) ||
        throw(ArgumentError("sigma_n_px must contain finite positive values."))

    has_length = haskey(data, "corr_length_px")
    has_matrix = haskey(data, "corr_matrix")
    has_length && has_matrix &&
        throw(ArgumentError("Provide at most one of corr_length_px and corr_matrix."))

    image_path = isabspath(image_file) ? image_file : joinpath(image_root, image_file)
    isfile(image_path) ||
        throw(ArgumentError("Edge image not found: $image_path"))
    isfile(expected_image_file) ||
        throw(ArgumentError("ArUco image not found: $expected_image_file"))
    expected_path = realpath(expected_image_file)
    resolved_image_path = realpath(image_path)
    image = FileIO.load(resolved_image_path)
    aruco_image = FileIO.load(expected_path)
    height, width = size(image, 1), size(image, 2)
    (height, width) == (size(aruco_image, 1), size(aruco_image, 2)) ||
        throw(ArgumentError("Edge image dimensions must match the ArUco image dimensions."))

    if has_length
        correlation_length = data["corr_length_px"]
        correlation_length isa Float64 && isfinite(correlation_length) &&
            correlation_length > 0 ||
            throw(ArgumentError("corr_length_px must be a finite positive Float64."))
        arc_length = zeros(Float64, n)
        for i in 2:n
            arc_length[i] = arc_length[i-1] + hypot(u[i] - u[i-1], v[i] - v[i-1])
        end
        correlation = [
            exp(-abs(arc_length[i] - arc_length[j]) / correlation_length)
            for i in 1:n, j in 1:n
        ]
        correlation_source = "corr_length_px"
    elseif has_matrix
        correlation = data["corr_matrix"]
        correlation isa Matrix{Float64} ||
            throw(ArgumentError("corr_matrix must be Matrix{Float64}."))
        correlation_source = "corr_matrix"
    else
        correlation = Matrix{Float64}(I, n, n)
        correlation_source = "uncorrelated"
    end

    correlation = _validate_correlation_matrix(correlation, n)
    all(x -> 1 <= x <= width, u) && all(y -> 1 <= y <= height, v) ||
        throw(ArgumentError("Edge points must lie inside the 1-based image bounds."))

    return (
        u=u,
        v=v,
        sigma_n_px=sigma_n_px,
        sigma_description=sigma_description,
        image_file=image_file,
        image_path=resolved_image_path,
        correlation=correlation,
        correlation_source=correlation_source,
    )
end


"""
Detect the first ArUco marker in an image without estimating its pose.
input:
- img_file [image path]
- K [3×3 camera matrix]
- dist [distortion coefficients]
- marker_dict [OpenCV ArUco dictionary or predefined dictionary ID]
- aruco_params [OpenCV ArUco detector parameters]
outputs:
- corners_px [4×2 matrix of detected image corners, px]
- ids [vector containing the detected marker ID]
"""
function detect_aruco_corners(img_file::String, K, dist; marker_dict=nothing, aruco_params=nothing)
    
    img_path = abspath(img_file)
    
    py"""
import cv2, numpy as np, os

img_path = $img_path
if not os.path.isfile(img_path):
    raise FileNotFoundError(f"Image file not found: {img_path}")
img = cv2.imread(img_path)
if img is None:
    raise ValueError(f"Could not read image: {img_path}")
if img.ndim == 2:
    gray = img
elif img.shape[2] == 4:
    gray = cv2.cvtColor(img, cv2.COLOR_RGBA2GRAY)
else:
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

aruco = cv2.aruco
selected_dict = $marker_dict
if selected_dict is None:
    selected_dict = aruco.getPredefinedDictionary(aruco.DICT_4X4_50)
elif isinstance(selected_dict, (int, np.integer)):
    selected_dict = aruco.getPredefinedDictionary(int(selected_dict))

detector_params = $aruco_params
if detector_params is None:
    detector_params = aruco.DetectorParameters()
    detector_params.adaptiveThreshWinSizeMin = 15
    detector_params.adaptiveThreshWinSizeMax = 99
    detector_params.adaptiveThreshWinSizeStep = 10
    detector_params.minMarkerPerimeterRate = 0.05
detector = aruco.ArucoDetector(selected_dict, detector_params)
detected_corners, detected_ids, _ = detector.detectMarkers(gray)
if detected_ids is None or len(detected_ids) == 0:
    raise ValueError(f"No ArUco marker detected in {img_path}")

corners_px = np.asarray(detected_corners[0], dtype=np.float64).reshape(4, 2)
ids = [int(np.asarray(detected_ids).reshape(-1)[0])]
result = (corners_px.tolist(), ids)
"""
    corners_raw, ids = py"result"
    corners_px = corners_raw isa AbstractMatrix ?
        Matrix{Float64}(corners_raw) :
        reduce(vcat, permutedims.(corners_raw))
    return Float64.(corners_px), Int.(ids)
end


function plot_aruco_corner_corrections(corner_results)

    """
    plot_aruco_corner_corrections(corner_results)

    Plot the x/y corrections returned by `fit_corner`, with their mean and
    68% and 95% scatter ellipses.

    Input:
    - `corner_results`: Results from `fit_corner`; 
                        `nothing` entries are ignored.
                        At least two valid results are required.

    Output:
    - A plot object showing the corner corrections in pixels.

    The ellipses summarise the scatter among valid corner corrections; they are
    not the individual covariance matrices returned by `fit_corner`.
    """

    valid_results = [result for result in corner_results if result !== nothing]
    length(valid_results) >= 2 ||
        throw(ArgumentError("At least two valid corner fits are needed to plot corrections."))

    corrections_px = Matrix{Float64}(undef, length(valid_results), 2)
    for (i, result) in enumerate(valid_results)
        corrections_px[i, :] = result.q
    end

    mean_correction_px = vec(mean(corrections_px; dims=1))
    centered_corrections_px = corrections_px .- mean_correction_px'
    correction_covariance_px2 =
        centered_corrections_px' * centered_corrections_px /
        (size(corrections_px, 1) - 1)

    covariance_eigen = eigen(Symmetric(correction_covariance_px2))
    angles = range(0, 2π; length=181)
    unit_circle = [cos.(angles)'; sin.(angles)']

    make_ellipse(scale) =
        reshape(mean_correction_px, 2, 1) .+
        scale .* covariance_eigen.vectors *
        Diagonal(sqrt.(max.(covariance_eigen.values, 0))) *
        unit_circle

    ellipse_68 = make_ellipse(sqrt(2.30))
    ellipse_95 = make_ellipse(sqrt(5.99))

    correction_plot = Plots.scatter(
        corrections_px[:, 1],
        corrections_px[:, 2];
        mc=:steelblue,
        ms=4,
        alpha=0.7,
        label="corner corrections",
        xlabel="x correction (px)",
        ylabel="y correction (px)",
        title="Edge-based correction relative to OpenCV",
        aspect_ratio=:equal,
        yflip=true,
    )
    Plots.vline!(correction_plot, [0.0]; lc=:black, ls=:dash, label="zero correction")
    Plots.hline!(correction_plot, [0.0]; lc=:black, ls=:dash, label="")
    Plots.plot!(
        correction_plot, ellipse_95[1, :], ellipse_95[2, :];
        lc=:orange, lw=2, label="95% scatter",
    )
    Plots.plot!(
        correction_plot, ellipse_68[1, :], ellipse_68[2, :];
        lc=:red, lw=2, label="68% scatter",
    )
    Plots.scatter!(
        correction_plot,
        [mean_correction_px[1]], [mean_correction_px[2]];
        mc=:red, ms=7, marker=:star5, label="mean correction",
    )

    return correction_plot
end


function plot_aruco_corner_diagnostic(intensity::AbstractMatrix{<:Real},
                                        corner_px::AbstractVector{<:Real},
                                        result;
                                        corner_index::Int,
                                        edge_index::Int=2,
                                        profile_index::Int=6,
                                    )

    """
    Create a side-by-side diagnostic plot of a fitted sigmoid intensity profile and
    its location on the image.

    Inputs:
    - `intensity`: 2D grayscale intensity matrix used for fitting.
    - `corner_px`: Corner coordinates `[x, y]` in the image's one-based pixel
    coordinate system.
    - `result`: Successful `fit_corner` result for this corner.
    - `corner_index`: Corner number, used in plot titles.
    - `edge_index`: Which fitted edge to show; defaults to `2`.
    - `profile_index`: Which profile on that edge to show; defaults to `6` and is
    clamped to the available profiles.

    Output:
    - A combined plot showing the sampled profile, sigmoid fit, transition region,
    and fitted corner overlaid on the image.

    The shaded transition region represents the sigmoid's intensity transition,
    not a confidence interval for the corner position.
    """

    1 <= edge_index <= length(result.edges) ||
        throw(ArgumentError("edge_index is out of range."))

    edge = result.edges[edge_index]
    isempty(edge.profiles) &&
        throw(ArgumentError("The selected edge has no usable profiles."))

    selected_profile_index = clamp(profile_index, 1, length(edge.profiles))
    profile = edge.profiles[selected_profile_index]
    parameters = profile.parameters

    edge_offset_px = parameters[3]
    transition_half_width_px = abs(parameters[4]) * log(39)
    transition_start_px = edge_offset_px - transition_half_width_px
    transition_end_px = edge_offset_px + transition_half_width_px

    fitted_d = range(minimum(result.d), maximum(result.d); length=400)
    sigmoid_plot = Plots.scatter(
        result.d,
        profile.profile_intensities;
        markersize=3,
        mc=:black,
        label="sampled intensity",
        xlabel="Offset across edge, d [px]",
        ylabel="Intensity",
        ylims=(0, 1),
        title="Corner $corner_index: edge $edge_index, profile $selected_profile_index",
    )
    Plots.plot!(
        sigmoid_plot, fitted_d, sigmoid(fitted_d, parameters);
        lc=:red, lw=2, label="sigmoid fit",
    )
    Plots.vspan!(
        sigmoid_plot, [transition_start_px, transition_end_px];
        fillcolor=:orange, fillalpha=0.25, label="95% transition interval",
    )
    Plots.vline!(
        sigmoid_plot, [edge_offset_px];
        lc=:blue, ls=:dash, label="chosen edge position",
    )
    Plots.scatter!(
        sigmoid_plot,
        [edge_offset_px], [sigmoid([edge_offset_px], parameters)[1]];
        mc=:blue, marker=:star5, label="sigmoid center",
    )

    profile_center = corner_px .+ profile.distance_along_edge .* edge.t
    edge_point = profile_center .+ edge_offset_px .* edge.n
    transition_start = edge_point .- transition_half_width_px .* edge.n
    transition_end = edge_point .+ transition_half_width_px .* edge.n
    fitted_corner = corner_px .+ result.q

    profile_positions = [
        item.distance_along_edge
        for fitted_edge in result.edges
        for item in fitted_edge.profiles
    ]
    plot_extent = ceil(Int, maximum(abs.(vcat(result.d, profile_positions))) + 4)

    image_height, image_width = size(intensity)
    center_x, center_y = round.(Int, corner_px)
    x_min = clamp(center_x - plot_extent, 1, image_width)
    x_max = clamp(center_x + plot_extent, 1, image_width)
    y_min = clamp(center_y - plot_extent, 1, image_height)
    y_max = clamp(center_y + plot_extent, 1, image_height)
    x_range = x_min:x_max
    y_range = y_min:y_max

    image_plot = Plots.heatmap(
        collect(x_range), collect(y_range), intensity[y_range, x_range];
        color=Plots.cgrad(:grays, rev=true),
        yflip=true,
        aspect_ratio=:equal,
        clims=(0, 1),
        colorbar=false,
        xlabel="x [px]",
        ylabel="y [px]",
        title="Corner $corner_index: fitted location",
    )

    sample_x = profile_center[1] .+ result.d .* edge.n[1]
    sample_y = profile_center[2] .+ result.d .* edge.n[2]
    Plots.plot!(image_plot, sample_x, sample_y; lc=:cyan, lw=1.5, label="sampled profile")
    Plots.plot!(
        image_plot,
        [transition_start[1], transition_end[1]],
        [transition_start[2], transition_end[2]];
        lc=:orange, lw=3, label="95% transition interval",
    )
    Plots.scatter!(
        image_plot, [edge_point[1]], [edge_point[2]];
        mc=:blue, marker=:star5, ms=6, label="sigmoid edge point",
    )
    Plots.scatter!(
        image_plot, [corner_px[1]], [corner_px[2]];
        mc=:magenta, ms=5, label="OpenCV corner",
    )
    Plots.scatter!(
        image_plot, [fitted_corner[1]], [fitted_corner[2]];
        mc=:lime, marker=:square, ms=4, label="fitted corner",
    )

    return Plots.plot(sigmoid_plot, image_plot; layout=(1, 2), size=(1100, 450))
end

function plot_aruco_edge_profile_residuals(
    intensity::AbstractMatrix{<:Real},
    result;
    corner_index::Int,
    edge_index::Int,
)
    1 <= edge_index <= length(result.edges) ||
        throw(ArgumentError("edge_index is out of range."))

    edge = result.edges[edge_index]
    corner_px = result.opencv_corner
    profiles = edge.profiles
    n_profiles = length(profiles)
    along_edge = [p.distance_along_edge for p in profiles]
    detected_offsets = [p.parameters[3] for p in profiles]

    # --- residuals ---
    norm_factor     = sqrt(1 + edge.b^2)                        # normalisation factor for perpendicular distance
    signed_distance = edge.residuals ./ norm_factor             # perpendicular distance from each detected edge point to the fitted line
    mean_distance = mean(signed_distance)                       # mean perpendicular distance from detected edge points to the fitted line
    σ_distance    = sqrt(edge.residual_variance) / norm_factor  # standard deviation of perpendicular distance from detected edge points to the fitted line

    # --- image-space geometry ---
    # unit normal to the fitted line, expressed in (s, c) coordinates
    line_normal_sc = [-edge.b, 1.0] ./ norm_factor
    to_pixels(s, c) = corner_px .+ s .* edge.t .+ c .* edge.n

    detected_points = [to_pixels(s, c) for (s, c) in zip(along_edge, detected_offsets)]
    foot_points = [
        to_pixels(
            s - d * line_normal_sc[1],
            c - d * line_normal_sc[2],
        ) for (s, c, d) in zip(along_edge, detected_offsets, signed_distance)
    ]

    line_positions = range(
        min(0.0, minimum(along_edge)),
        max(0.0, maximum(along_edge));
        length=100,
    )
    fitted_line = [to_pixels(s, edge.a0 + edge.b * s) for s in line_positions]

    # --- image panel ---
    plot_extent = ceil(Int, maximum(abs.(vcat(along_edge, result.d))) + 8)
    image_height, image_width = size(intensity)
    center_x, center_y = round.(Int, corner_px)
    x_range = clamp(center_x - plot_extent, 1, image_width):clamp(center_x + plot_extent, 1, image_width)
    y_range = clamp(center_y - plot_extent, 1, image_height):clamp(center_y + plot_extent, 1, image_height)

    image_plot = Plots.heatmap(
        collect(x_range), collect(y_range), intensity[y_range, x_range];
        color=Plots.cgrad(:grays, rev=true),
        yflip=true, aspect_ratio=:equal, clims=(0, 1), colorbar=false,
        xlabel="x [px]", ylabel="y [px]",
        title="Corner $corner_index, edge $edge_index",
    )
    Plots.plot!(image_plot,
        [p[1] for p in fitted_line], [p[2] for p in fitted_line];
        lc=:lime, lw=2, label="fitted edge")
    for i in 1:n_profiles
        Plots.plot!(image_plot,
            [foot_points[i][1], detected_points[i][1]],
            [foot_points[i][2], detected_points[i][2]];
            lc=:orange, lw=2,
            label=i == 1 ? "perpendicular distance" : "")
    end
    Plots.scatter!(image_plot,
        [p[1] for p in detected_points], [p[2] for p in detected_points];
        mc=:steelblue, ms=4, label="detected edge points")

    # --- distance panel ---
    y_half = 1.3 * max(maximum(abs.(signed_distance .- mean_distance)), σ_distance)
    y_limits = (mean_distance - y_half, mean_distance + y_half)

    profile_numbers = collect(1:n_profiles)
    distance_plot = Plots.hspan(
        [mean_distance - σ_distance, mean_distance + σ_distance];
        fillcolor=:orange, fillalpha=0.25, linecolor=:transparent,
        label="mean ± 1σ (σ = $(round(σ_distance; digits=3)) px)",
        xlims=(0.5, n_profiles + 0.5), ylims=y_limits,
        xticks=(profile_numbers, ["P$i" for i in profile_numbers]),
        xlabel="Profile",
        ylabel="Perpendicular distance to fitted edge [px]",
        title="Edge scatter about fitted line",
    )
    Plots.hline!(distance_plot, [mean_distance];
        lc=:black, lw=2, label="mean = $(round(mean_distance; digits=3)) px")
    Plots.hline!(distance_plot, [mean_distance - σ_distance, mean_distance + σ_distance];
        lc=:orange, ls=:dash, lw=1.5, label="")
    Plots.scatter!(distance_plot, profile_numbers, signed_distance;
        mc=:steelblue, ms=6, label="detected edge")

    return Plots.plot(image_plot, distance_plot;
        layout=(1, 2), size=(1200, 500))
end


# Image selection

"""
List image files in a directory whose names match the requested extensions.

input:
- img_dir [directory path]
- image_extensions [file extensions]

outputs:
- Sorted absolute image paths [paths]
"""
function list_camera_images(
    img_dir::String;
    image_extensions::Vector{String} = String[".JPG"],
)::Vector{String}
    isdir(img_dir) || return String[]
    files = String[]
    for name in readdir(img_dir)
        path = joinpath(img_dir, name)
        if isfile(path) && any(endswith(lowercase(name), lowercase(ext)) for ext in image_extensions)
            push!(files, abspath(path))
        end
    end
    return sort(files)
end

"""
Choose a random subset of calibration images.

input:
- img_dir [directory path]
- subset_size [images]
- seed [integer or `nothing`]
- image_extensions [file extensions]

outputs:
- Selected absolute image paths [paths]
"""
function random_image_subset(
    img_dir::String;
    subset_size::Int                  = 5,
    seed::Union{Int,Nothing}          = nothing,
    image_extensions::Vector{String}  = String[".JPG"],
)::Vector{String}
    files = list_camera_images(img_dir; image_extensions)
    (isempty(files) || subset_size <= 0) && return String[]
    shuffled = copy(files)
    rng = seed !== nothing ? MersenneTwister(seed) : Random.GLOBAL_RNG
    Random.shuffle!(rng, shuffled)
    return shuffled[1:min(subset_size, length(shuffled))]
end

# Calibration

"""
Calibrate camera intrinsics from checkerboard images using OpenCV.

input:
- image_files [image paths]
- checkerboard_size [inner corners]
- square_size [m]
- out_dir [directory path or `nothing`]
- reference_K [camera matrix or `nothing`]
- reference_dist [distortion coefficients or `nothing`]

outputs:
- Calibration results Dict [camera matrix, distortion, RMS in px, image paths, count, success]
- std_intrinsics [σ per parameter: fx fy cx cy k1 k2 p1 p2 k3, px of dimensieloos]
"""
function calibrate_camera_from_images(
    image_files::Vector{String};
    checkerboard_size::Tuple{Int,Int}   = (13, 9),
    square_size::Float64                = 0.021,
    out_dir::Union{String,Nothing}      = nothing,
    reference_K                         = nothing,
    reference_dist                      = nothing,
)
    """
    Build the fallback result used when calibration cannot run.

    input:
    - None

    outputs:
    - Failure result Dict [empty calibration result]
    """
    _fail() = Dict("camera_matrix" => nothing, "dist_coeffs" => nothing,
                   "rms" => nothing, "images" => String[], "n_images" => 0,
                   "success" => false, "std_intrinsics" => Float64[])

    isempty(image_files) && return _fail()
    out_dir !== nothing && ensure_directory(out_dir)

    result = nothing
    try
        result = py"""
        import cv2, numpy as np, os

        checker     = tuple($checkerboard_size)
        square      = float($square_size)
        image_files = $image_files
        out_dir     = None if $out_dir is None else $out_dir

        if out_dir is not None:
            os.makedirs(out_dir, exist_ok=True)

        objpoints, imgpoints, used_images = [], [], []

        objp = np.zeros((checker[0]*checker[1], 3), np.float32)
        objp[:, :2] = np.mgrid[0:checker[0], 0:checker[1]].T.reshape(-1, 2)
        objp *= square

        last_gray = None

        for img_path in image_files:
            img = cv2.imread(img_path)
            if img is None:
                continue
            if img.ndim == 2:
                img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
            elif img.shape[2] == 4:
                img = cv2.cvtColor(img, cv2.COLOR_RGBA2BGR)

            gray      = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
            last_gray = gray

            ret, corners = cv2.findChessboardCornersSB(gray, checker)
            if not ret:
                continue

            corners = corners.astype(np.float32)
            criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 30, 0.001)
            corners = cv2.cornerSubPix(gray, corners, (11, 11), (-1, -1), criteria)

            objpoints.append(objp)
            imgpoints.append(corners)
            used_images.append(os.path.basename(img_path))

            if out_dir is not None:
                cv2.drawChessboardCorners(img, checker, corners, True)
                cv2.imwrite(os.path.join(out_dir, os.path.basename(img_path)), img)

        if len(objpoints) > 0 and last_gray is not None:
            h, w    = last_gray.shape[:2]
            init_K  = None if $reference_K   is None else np.array($reference_K,   dtype=np.float64)
            init_d  = None if $reference_dist is None else np.array($reference_dist, dtype=np.float64).reshape(1, -1)
            flags   = cv2.CALIB_USE_INTRINSIC_GUESS if init_K is not None else 0
            ret, mtx, dist, rvecs, tvecs, std_int, std_ext, per_view_errors = \
                cv2.calibrateCameraExtended(objpoints, imgpoints, (w, h), init_K, init_d, flags=flags)
            result = (float(ret), mtx.tolist(), dist.reshape(-1).tolist(), used_images,
                      len(used_images), True, [o.tolist() for o in objpoints],
                      [i.reshape(-1,2).tolist() for i in imgpoints],
                      [r.reshape(-1).tolist() for r in rvecs],
                      [t.reshape(-1).tolist() for t in tvecs], (w, h),
                      std_int.reshape(-1).tolist())
        else:
            result = (None, None, None, used_images, len(used_images), False)

        result
        """
        result = py"result"
    catch err
        @warn "Calibration failed: $err"
        return _fail()
    end

    result === nothing && return _fail()
    !(result isa Tuple && length(result) == 12) && return _fail()

    rms, cammat, distcoef, used_images, n_images, success,
    objpoints, imgpoints, rvecs, tvecs, image_size, std_intrinsics = result
    return Dict("camera_matrix" => cammat, "dist_coeffs" => distcoef,
                "rms" => rms, "images" => used_images,
                "n_images" => n_images, "success" => success,
                "objpoints" => objpoints, "imgpoints" => imgpoints,
                "rvecs" => rvecs, "tvecs" => tvecs, "image_size" => image_size,
                "std_intrinsics" => Float64.(std_intrinsics[1:9]))
end


"""
Plot pooled x/y reprojection residuals and save the scatter plot.

input:
- residuals [N×2 array, px]
- save_dir [output directory path]

outputs:
- p [scatter plot]
- PNG file [saved in `save_dir`]
"""
function plot_reprojection_residuals(residuals; save_dir::String = ".")
    size(residuals, 2) == 2 || throw(ArgumentError("residuals must have two columns (x and y)."))
    mkpath(save_dir)
    p = scatter(residuals[:, 1], residuals[:, 2]; marker=:circle, markersize=2,
                alpha=0.45, xlabel="x residual [px]", ylabel="y residual [px]",
                title="Pooled checkerboard reprojection residuals", legend=false)
    savefig(p, joinpath(save_dir, "reprojection_residuals.png"))
    return p
end


"""
Plot intrinsic-parameter subset statistics against the all-image reference.

input:
- summary [output from `summarise_by_subset`]
- K_ref [reference camera matrix]
- dist_ref [reference distortion coefficients]
- save_dir [output directory path]

outputs:
- p_intrinsic [intrinsic-parameter plot]
- p_distortion [distortion-coefficient plot]
- PNG files [saved in `save_dir`]
"""
function plot_subset_influence(
    summary ::Dict{String,NamedTuple},
    K_ref   ::Matrix{Float64},
    dist_ref::Vector{Float64};
    save_dir::String = ".",
)
    if !any(s -> !isempty(s.n), values(summary))
        @warn "No successful subset calibrations; skipping subset influence plots."
        return nothing, nothing
    end
    if size(K_ref) != (3, 3) || length(dist_ref) < 5 ||
       !all(isfinite, K_ref) || !all(isfinite, dist_ref[1:5])
        @warn "Invalid reference calibration; skipping subset influence plots."
        return nothing, nothing
    end

    mkpath(save_dir)
    ref_vals = Dict("fx" => K_ref[1,1], "fy" => K_ref[2,2],
                    "cx" => K_ref[1,3], "cy" => K_ref[2,3],
                    "k1" => dist_ref[1], "k2" => dist_ref[2],
                    "p1" => dist_ref[3], "p2" => dist_ref[4], "k3" => dist_ref[5])

    """
    Build a parameter plot with its reference value.

    input:
    - name [intrinsic parameter name]

    outputs:
    - p [parameter plot]
    """
    function _panel(name)
        s = summary[name]
        p = plot(s.n, s.mean; ribbon=s.std, fillalpha=0.25, marker=:circle,
                 xlabel="# calibration images", ylabel="value", title=name, legend=false)
        hline!(p, [ref_vals[name]]; linestyle=:dash, color=:black)
        return p
    end

    p_int  = plot([_panel(n) for n in ["fx","fy","cx","cy"]]...,
                  layout=(2,2), size=(900,700),
                  plot_title="Intrinsic parameters vs. # calibration images")
    p_dist = plot([_panel(n) for n in ["k1","k2","p1","p2","k3"]]...,
                  layout=(2,3), size=(1100,700),
                  plot_title="Distortion coefficients vs. # calibration images")

    savefig(p_int,  joinpath(save_dir, "subset_intrinsics.png"))
    savefig(p_dist, joinpath(save_dir, "subset_distortion.png"))

    return p_int, p_dist
end


"""
Plot each successful calibration sample against subset size, with color identifying its seed.

input:
- results [subset-study results Dict]
- K_ref [reference camera matrix]
- dist_ref [reference distortion coefficients]
- save_dir [output directory path]

outputs:
- p_intrinsic [per-seed intrinsic-parameter plot]
- p_distortion [per-seed distortion-coefficient plot]
- PNG files [saved in `save_dir`]
"""
function plot_subset_seed_values(
    results::Dict{Int,Dict},
    K_ref::Matrix{Float64},
    dist_ref::Vector{Float64};
    save_dir::String = ".",
)
    isempty(results) && return nothing, nothing
    if size(K_ref) != (3, 3) || length(dist_ref) < 5 ||
       !all(isfinite, K_ref) || !all(isfinite, dist_ref[1:5])
        @warn "Invalid reference calibration; skipping per-seed plots."
        return nothing, nothing
    end

    mkpath(save_dir)
    ref_vals = Dict("fx" => K_ref[1,1], "fy" => K_ref[2,2],
                    "cx" => K_ref[1,3], "cy" => K_ref[2,3],
                    "k1" => dist_ref[1], "k2" => dist_ref[2],
                    "p1" => dist_ref[3], "p2" => dist_ref[4], "k3" => dist_ref[5])
    accessors = Dict(
        "fx" => (K, d) -> K[1,1], "fy" => (K, d) -> K[2,2],
        "cx" => (K, d) -> K[1,3], "cy" => (K, d) -> K[2,3],
        "k1" => (K, d) -> d[1], "k2" => (K, d) -> d[2],
        "p1" => (K, d) -> d[3], "p2" => (K, d) -> d[4],
        "k3" => (K, d) -> d[5],
    )

    subset_sizes = sort(collect(keys(results)))
    seed_ids = sort(unique(vcat([
        get(results[n], "seed_ids", collect(1:length(results[n]["K_samples"])))
        for n in subset_sizes
    ]...)))
    seed_indices = Dict(seed => i for (i, seed) in enumerate(seed_ids))
    seed_clims = (0.5, length(seed_ids) + 0.5)

    function _seed_panel(name)
        x_values, y_values, color_values = Int[], Float64[], Int[]
        for seed in seed_ids
            for n in subset_sizes
                data = results[n]
                seeds = get(data, "seed_ids", collect(1:length(data["K_samples"])))
                sample_index = findfirst(==(seed), seeds)
                sample_index === nothing && continue
                K_samples = data["K_samples"]
                dist_samples = data["dist_samples"]
                sample_index <= length(K_samples) && sample_index <= length(dist_samples) || continue
                length(dist_samples[sample_index]) >= 5 || continue
                push!(x_values, n)
                push!(y_values, accessors[name](K_samples[sample_index], dist_samples[sample_index]))
                push!(color_values, seed_indices[seed])
            end
        end
        p = plot(; xlabel="# calibration images", ylabel="value",
                 title=name, legend=:best, xticks=subset_sizes)
        if !isempty(x_values)
            scatter!(p, x_values, y_values;
                     marker=:circle, marker_z=color_values, markercolor=:rainbow,
                     clims=seed_clims, colorbar=false, label="")
            scatter!(p, [NaN], [NaN]; marker=:circle, markercolor=:white,
                     markerstrokecolor=:black, markerstrokewidth=1.5,
                     colorbar=false, label="Seeds")
        end
        hline!(p, [ref_vals[name]]; linestyle=:dash, color=:black, label="all images")
        return p
    end

    intrinsic_names = ["fx", "fy", "cx", "cy"]
    distortion_names = ["k1", "k2", "p1", "p2", "k3"]
    p_intrinsic = plot([_seed_panel(name) for name in intrinsic_names]...;
                       layout=(2,2), size=(900,700),
                       plot_title="Intrinsic parameters by subset size and seed")
    p_distortion = plot([_seed_panel(name) for name in distortion_names]...;
                        layout=(2,3), size=(1100,700),
                        plot_title="Distortion coefficients by subset size and seed")

    savefig(p_intrinsic, joinpath(save_dir, "subset_intrinsics_by_seed.png"))
    savefig(p_distortion, joinpath(save_dir, "subset_distortion_by_seed.png"))
    return p_intrinsic, p_distortion
end


# Bootstrap plotting

"""
Plot histograms of bootstrap rotation and translation samples.
input:
- R_samples [rotation matrices, one 3×3 matrix per sample]
- t_samples [translation vectors, one 3-vector per sample, mm]
- save_fig [whether to save the plot]
- out_dir [output directory path]
- name_suffix [optional suffix for the saved plot filename]
outputs:
- plt [plot object]
"""
function plot_extrinsic_distributions(
    R_samples::Vector{Matrix{Float64}},
    t_samples::Vector{Vector{Float64}};
    save_fig::Bool = false,
    out_dir::String = ".",
    name_suffix::String = "",
)
    length(R_samples) == length(t_samples) ||
        throw(ArgumentError("R_samples and t_samples must have the same number of samples."))
    isempty(R_samples) &&
        throw(ArgumentError("At least one pose sample is required to plot distributions."))
    all(size(R) == (3, 3) && all(isfinite, R) for R in R_samples) ||
        throw(ArgumentError("Every rotation sample must be a finite 3×3 matrix."))
    all(length(t) == 3 && all(isfinite, t) for t in t_samples) ||
        throw(ArgumentError("Every translation sample must contain three finite values."))

    py"""
import cv2, numpy as np

rotations = $R_samples
rvec_samples = []
for rotation in rotations:
    rvec, _ = cv2.Rodrigues(np.array(rotation, dtype=np.float64))
    rvec_samples.append(rvec.reshape(-1).tolist())
mean_rvec = np.mean(np.asarray(rvec_samples, dtype=np.float64), axis=0).reshape(3, 1)
mean_rotation, _ = cv2.Rodrigues(mean_rvec)
delta_rvec_samples = []
for rotation in rotations:
    delta_rotation = mean_rotation.T @ np.array(rotation, dtype=np.float64)
    delta_rvec, _ = cv2.Rodrigues(delta_rotation)
    delta_rvec_samples.append(delta_rvec.reshape(-1).tolist())
result = (delta_rvec_samples, mean_rotation.tolist())
"""
    delta_rvec_raw, _ = py"result"
    delta_rvec_samples_deg = rad2deg.(Matrix{Float64}(delta_rvec_raw))
    translation_values_mm = reduce(vcat, permutedims.(t_samples))
    parameter_samples = (
        delta_rvec_samples_deg[:, 1], delta_rvec_samples_deg[:, 2],
        delta_rvec_samples_deg[:, 3],
        translation_values_mm[:, 1], translation_values_mm[:, 2], translation_values_mm[:, 3],
    )
    parameter_names = ("Δrx", "Δry", "Δrz", "tx", "ty", "tz")
    parameter_units = ("deg", "deg", "deg", "mm", "mm", "mm")

    plt = plot(layout=(2, 3), size=(1200, 700))
    for i in eachindex(parameter_names)
        values = parameter_samples[i]
        μ = mean(values)
        σ = std(values)
        histogram!(plt[i], values;
                   bins=30, title=parameter_names[i],
                   xlabel="[$(parameter_units[i])]", ylabel="Count",
                   legend=false)
        vline!(plt[i], [μ]; color=:red, label="")
        vline!(plt[i], [μ - σ, μ + σ]; color=:gray, linestyle=:dash, label="")
    end

    if save_fig
        mkpath(out_dir)
        suffix = isempty(name_suffix) ? "" : "_$(name_suffix)"
        savefig(plt, joinpath(out_dir, "extrinsic_distributions$(suffix).png"))
    end
    return plt
end


"""
Plot stability of extrinsic bootstrap standard deviations.
input:
- R_samples [vector of N 3×3 rotation matrices]
- t_samples [vector of N 3-vectors, mm]
- checkpoints [Vector{Int}, number of samples per checkpoint]
- save_fig [Bool]
- out_dir [directory path]
- name_suffix [String appended to the filename]
outputs:
- plt [Plots object]
- PNG file [saved in out_dir as stability_extrinsics_<name_suffix>.png]
"""
function plot_extrinsic_bootstrap_stability(
    R_samples::Vector{Matrix{Float64}},
    t_samples::Vector{Vector{Float64}};
    checkpoints::Vector{Int},
    save_fig::Bool = false,
    out_dir::String = ".",
    name_suffix::String = "",
)
    length(R_samples) == length(t_samples) ||
        throw(ArgumentError("R_samples and t_samples must have the same number of samples."))
    n_samples = length(R_samples)
    n_samples >= 2 ||
        throw(ArgumentError("At least two pose samples are required to calculate standard deviations."))
    isempty(checkpoints) &&
        throw(ArgumentError("checkpoints must contain at least one sample count."))
    all(2 .<= checkpoints) && all(checkpoints .<= n_samples) ||
        throw(ArgumentError("Every checkpoint must be between 2 and the number of samples ($n_samples)."))
    issorted(checkpoints) && length(unique(checkpoints)) == length(checkpoints) ||
        throw(ArgumentError("checkpoints must be strictly increasing."))
    all(size(R) == (3, 3) && all(isfinite, R) for R in R_samples) ||
        throw(ArgumentError("Every rotation sample must be a finite 3×3 matrix."))
    all(length(t) == 3 && all(isfinite, t) for t in t_samples) ||
        throw(ArgumentError("Every translation sample must contain three finite values."))

    py"""
import cv2, numpy as np

rotations = $R_samples
rvec_samples = []
for rotation in rotations:
    rvec, _ = cv2.Rodrigues(np.array(rotation, dtype=np.float64))
    rvec_samples.append(rvec.reshape(-1).tolist())
rotvecs = np.asarray(rvec_samples, dtype=np.float64)
checkpoint_values = $checkpoints
std_rotation_per_checkpoint = []
for sample_count in checkpoint_values:
    mean_rvec = np.mean(rotvecs[:sample_count], axis=0).reshape(3, 1)
    mean_rotation, _ = cv2.Rodrigues(mean_rvec)
    delta_rvecs = []
    for rotation in rotations[:sample_count]:
        delta_rotation = mean_rotation.T @ np.array(rotation, dtype=np.float64)
        delta_rvec, _ = cv2.Rodrigues(delta_rotation)
        delta_rvecs.append(delta_rvec.reshape(-1).tolist())
    std_rotation_per_checkpoint.append(
        np.std(np.asarray(delta_rvecs, dtype=np.float64), axis=0, ddof=1).tolist())
result = std_rotation_per_checkpoint
"""
    rotation_std_per_checkpoint = rad2deg.(Matrix{Float64}(py"result"))
    t_matrix = reduce(vcat, permutedims.(t_samples))
    translation_std_per_checkpoint = [
        std(t_matrix[1:M, i]) for M in checkpoints, i in 1:3
    ]
    std_per_checkpoint = hcat(rotation_std_per_checkpoint, translation_std_per_checkpoint)
    std_final = vec(std_per_checkpoint[end, :])
    parameter_names = ("Δrx", "Δry", "Δrz", "tx", "ty", "tz")
    parameter_units = ("deg", "deg", "deg", "mm", "mm", "mm")

    plt = plot(layout=(2, 3), size=(1200, 700))
    for i in 1:6
        plot!(plt[i], checkpoints, std_per_checkpoint[:, i];
              marker=:circle, color=:blue, label="std",
              xlabel="N_boot samples",
              ylabel="std [$(parameter_units[i])]",
              title=parameter_names[i], legend=:best,
              xticks=checkpoints)
        hline!(plt[i], [std_final[i]];
               linestyle=:dash, color=:black, label="dashed")
    end

    if save_fig
        mkpath(out_dir)
        suffix = isempty(name_suffix) ? "" : "_$(name_suffix)"
        savefig(plt, joinpath(out_dir, "stability_extrinsics$(suffix).png"))
    end
    return plt
end


"""
Plot convergence of the mean of extrinsic bootstrap parameters.
input:
- R_samples [vector of N 3×3 rotation matrices]
- t_samples [vector of N 3-vectors, mm]
- checkpoints [Vector{Int}]
- n_seeds [number of permutations]
- save_fig [Bool]
- out_dir [directory path]
- name_suffix [String appended to the filename]
outputs:
- plt [Plots object]
- PNG file [saved in out_dir as convergence_mean_extrinsics_<name_suffix>.png]
"""
function plot_extrinsic_convergence_mean(
    R_samples::Vector{Matrix{Float64}},
    t_samples::Vector{Vector{Float64}};
    checkpoints::Vector{Int},
    n_seeds::Int,
    save_fig::Bool = false,
    out_dir::String = ".",
    name_suffix::String = "",
)
    length(R_samples) == length(t_samples) ||
        throw(ArgumentError("R_samples and t_samples must have the same number of samples."))
    n_samples = length(R_samples)
    n_samples >= 2 ||
        throw(ArgumentError("At least two pose samples are required to estimate convergence."))
    n_seeds >= 2 ||
        throw(ArgumentError("n_seeds must be at least 2 to estimate spread across permutations."))
    isempty(checkpoints) &&
        throw(ArgumentError("checkpoints must contain at least one sample count."))
    all(1 .<= checkpoints) && all(checkpoints .<= n_samples) ||
        throw(ArgumentError("Every checkpoint must be between 1 and the number of samples ($n_samples)."))
    issorted(checkpoints) && length(unique(checkpoints)) == length(checkpoints) ||
        throw(ArgumentError("checkpoints must be strictly increasing."))
    all(size(R) == (3, 3) && all(isfinite, R) for R in R_samples) ||
        throw(ArgumentError("Every rotation sample must be a finite 3×3 matrix."))
    all(length(t) == 3 && all(isfinite, t) for t in t_samples) ||
        throw(ArgumentError("Every translation sample must contain three finite values."))

    py"""
import cv2, numpy as np

rotations = $R_samples
rvec_samples = []
for rotation in rotations:
    rvec, _ = cv2.Rodrigues(np.array(rotation, dtype=np.float64))
    rvec_samples.append(rvec.reshape(-1).tolist())
mean_rvec = np.mean(np.asarray(rvec_samples, dtype=np.float64), axis=0).reshape(3, 1)
mean_rotation, _ = cv2.Rodrigues(mean_rvec)
delta_rvec_samples = []
for rotation in rotations:
    delta_rotation = mean_rotation.T @ np.array(rotation, dtype=np.float64)
    delta_rvec, _ = cv2.Rodrigues(delta_rotation)
    delta_rvec_samples.append(delta_rvec.reshape(-1).tolist())
result = delta_rvec_samples
"""
    rotation_values_deg = rad2deg.(Matrix{Float64}(py"result"))
    translation_values_mm = reduce(vcat, permutedims.(t_samples))
    parameter_samples = hcat(rotation_values_deg, translation_values_mm)
    n_checkpoints = length(checkpoints)
    mean_per_checkpoint = zeros(Float64, n_checkpoints, 6)
    spread_per_checkpoint = zeros(Float64, n_checkpoints, 6)
    shuffled_indices = [
        Random.shuffle(Random.MersenneTwister(seed), collect(1:n_samples))
        for seed in 1:n_seeds
    ]

    for (checkpoint_index, sample_count) in enumerate(checkpoints)
        seed_means = Matrix{Float64}(undef, n_seeds, 6)
        for seed_index in 1:n_seeds
            indices = view(shuffled_indices[seed_index], 1:sample_count)
            seed_means[seed_index, :] .=
                vec(mean(parameter_samples[indices, :]; dims=1))
        end
        mean_per_checkpoint[checkpoint_index, :] .=
            vec(mean(seed_means; dims=1))
        spread_per_checkpoint[checkpoint_index, :] .=
            vec(std(seed_means; dims=1))
    end

    mean_final = vec(mean(parameter_samples; dims=1))
    parameter_names = ("Δrx", "Δry", "Δrz", "tx", "ty", "tz")
    parameter_units = ("deg", "deg", "deg", "mm", "mm", "mm")
    latex_strings = Base.require(@__MODULE__, :LaTeXStrings)
    x_min = checkpoints[1]
    x_max = checkpoints[end]
    tick_exponents = 1:floor(Int, log10(x_max))
    tick_positions = [10^i for i in tick_exponents]
    tick_labels = Any[
        Base.invokelatest(latex_strings.LaTeXString, "10^$(i)")
        for i in tick_exponents
    ]
    if x_max ∉ tick_positions
        push!(tick_positions, x_max)
        push!(tick_labels, string(x_max))
    end
    xlims = (x_min, x_max * 1.05)
    xticks = (tick_positions, tick_labels)
    plt = plot(layout=(2, 3), size=(1200, 700))
    for i in 1:6
        plot!(plt[i], checkpoints, mean_per_checkpoint[:, i];
              ribbon=spread_per_checkpoint[:, i], marker=:circle,
              xlabel="N_boot samples", ylabel="mean [$(parameter_units[i])]",
              title=parameter_names[i], label="mean ± std over permutations",
              xticks=xticks, xlims=xlims,
              xaxis=:log, yaxis=:identity)
        C = spread_per_checkpoint[1, i] * sqrt(checkpoints[1])
        reference_offset = C ./ sqrt.(checkpoints)
        plot!(plt[i], checkpoints, mean_final[i] .+ reference_offset;
              color=:red, linestyle=:dash, label="mean_final ± C/√N")
        plot!(plt[i], checkpoints, mean_final[i] .- reference_offset;
              color=:red, linestyle=:dash, label="")
    end

    if save_fig
        mkpath(out_dir)
        suffix = isempty(name_suffix) ? "" : "_$(name_suffix)"
        Base.invokelatest(
            savefig, plt, joinpath(out_dir, "convergence_mean_extrinsics$(suffix).png"),
        )
    end
    return plt
end


"""
Save extrinsic pose means and covariance matrices to a JLD2 file.
input:
- result [tuple, NamedTuple, or Dict from extrinsic_covariance]
- output_path [JLD2 file path]
outputs:
- success [whether the file was saved]
"""
function save_extrinsic_uncertainty(result, output_path::String)
    if result isa Tuple && length(result) == 4
        Σ_t, Σ_rvec, rvec_mean, t_mean = result
    elseif result isa NamedTuple
        all(name -> haskey(result, name), (:Σ_t, :Σ_rvec, :rvec_mean, :t_mean)) ||
            throw(ArgumentError("result must contain Σ_t, Σ_rvec, rvec_mean, and t_mean."))
        Σ_t, Σ_rvec, rvec_mean, t_mean =
            result.Σ_t, result.Σ_rvec, result.rvec_mean, result.t_mean
    elseif result isa AbstractDict
        get_result_value(name) = haskey(result, name) ? result[name] :
            haskey(result, String(name)) ? result[String(name)] :
            throw(ArgumentError("result must contain $name."))
        Σ_t, Σ_rvec, rvec_mean, t_mean =
            get_result_value(:Σ_t), get_result_value(:Σ_rvec),
            get_result_value(:rvec_mean), get_result_value(:t_mean)
    else
        throw(ArgumentError("result must be a four-value tuple, NamedTuple, or Dict."))
    end

    Σ_t = Matrix{Float64}(Σ_t)
    Σ_rvec = Matrix{Float64}(Σ_rvec)
    rvec_mean = Vector{Float64}(rvec_mean)
    t_mean = Vector{Float64}(t_mean)
    size(Σ_t) == (3, 3) && size(Σ_rvec) == (3, 3) ||
        throw(ArgumentError("Σ_t and Σ_rvec must both be 3×3 matrices."))
    length(rvec_mean) == 3 && length(t_mean) == 3 ||
        throw(ArgumentError("rvec_mean and t_mean must each contain three values."))
    all(isfinite, Σ_t) && all(isfinite, Σ_rvec) &&
        all(isfinite, rvec_mean) && all(isfinite, t_mean) ||
        throw(ArgumentError("Extrinsic uncertainty results must contain only finite values."))
    isempty(output_path) && throw(ArgumentError("output_path must not be empty."))

    try
        output_dir = dirname(abspath(output_path))
        mkpath(output_dir)
        jldsave(output_path; Σ_t, Σ_rvec, rvec_mean, t_mean)
        return true
    catch err
        @error "Failed to save extrinsic uncertainty to $output_path" exception=(err, catch_backtrace())
        return false
    end
end

function make_normal_band(r, z, normals, σ)
    upper = [Float64[r[i] + σ[i] * normals[i, 1],
                     z[i] + σ[i] * normals[i, 2]] for i in eachindex(r)]
    lower = [Float64[r[i] - σ[i] * normals[i, 1],
                     z[i] - σ[i] * normals[i, 2]] for i in eachindex(r)]
    polygon = vcat(upper, reverse(lower))
    return Shape(first.(polygon), last.(polygon))
end
