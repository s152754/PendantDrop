project_root = normpath(joinpath(@__DIR__, ".."))
include(joinpath(project_root, "src", "camera.jl"))
using .Camera
using FileIO
using Images
using Interpolations
using JLD2
using LsqFit
using LinearAlgebra
using Plots
using Printf
using Statistics

# ── configuration ─────────────────────────────────────────────────────────────

zoom_level = "zoom5"        # Zoom setting of the droplet image 

# File paths for calibration and output
intrinsics_file = joinpath(
    project_root, "data", "calibration_camera", zoom_level, "detected", "camint_all.jld2",
)
out_dir = joinpath(project_root, "data", "test_droplets")
img_file = joinpath(out_dir, "$(zoom_level).JPG")

# ArUco marker length in millimeters (used for extrinsic calibration)
marker_length_mm = 10.0

# Settings for bootstrap sampling and convergence/stability analysis
N_boot = 1000               # Number of bootstrap samples for extrinsic uncertainty estimation
checkpoints_convergence = unique(round.(Int,
                                10 .^ range(log10(5), 
                                log10(N_boot), 
                                length=20)))            
checkpoints_stability = filter(
                                x -> x <= N_boot,
                                [100, 200, 400, 800, 1500, 2500, 3500, 5000],
                            )
n_seeds_convergence = 10    # Number of random seeds to test convergence of the mean extrinsic parameters
seed = 123                  # Random seed for reproducibility of bootstrap sampling

# Settings for edge profile sampling and fitting. These are tuned to the ArUco marker size and image resolution.
edge_distance_fraction = 0.4    # fraction of shortest edge to use for edge profiles (must be > profile_start_offset + min_profile_extent)
n_profiles = 6                  # number of edge profiles to sample per corner
profile_start_offset = 8.0      # min distance from corner to start of edge profile (in pixels)
min_profile_extent = 5.0        # min distance from start to end of edge profile (in pixels)
edge_sel = 2                    # which edge to plot for diagnostics (1-4, counter-clockwise from top-left)
profile_sel = 6                 # which profile to plot for diagnostics (1-n_profiles, from corner outwards)

# ── load intrinsics ────────────────────────────────────────────────────────────

isfile(intrinsics_file) || error("Intrinsic calibration file not found: $intrinsics_file")
intrinsics = JLD2.load(intrinsics_file)
K_raw = get(intrinsics, "camera_matrix", get(intrinsics, "K", nothing))
dist_raw = get(intrinsics, "dist_coeffs", get(intrinsics, "dist", nothing))
K_raw === nothing && error("No camera matrix found in $intrinsics_file")
dist_raw === nothing && error("No distortion coefficients found in $intrinsics_file")
K = Matrix{Float64}(K_raw)
dist = vec(Float64.(dist_raw))

# ── load or compute bootstrap samples ─────────────────────────────────────────

mkpath(out_dir)
isfile(img_file) || error("Droplet image not found: $img_file")
corners_px, marker_ids = Camera.IntrinsicUncertainty.detect_aruco_corners(
    img_file, K, dist,
)
println("Detected ArUco marker ID(s): ", join(marker_ids, ", "))
size(corners_px) == (4, 2) ||
    error("Expected four corners for one ArUco marker; got $(size(corners_px)).")

# ── compute intensity interpolant ─────────────────────────────────────────────
image = FileIO.load(img_file)
intensity = 1 .- Float64.(Gray.(image))
intensity_interpolant = linear_interpolation(
    (1:size(intensity, 1), 1:size(intensity, 2)),
    intensity;
    extrapolation_bc=NaN,
)

# ── fit edge profiles and compute corner uncertainties ─────────────────────────────

# OpenCV corners are zero-based; shift to Julia's one-based image coordinates while sampling/fitting
# Correction q is in pixels + applies unchanged to the original OpenCV coordinates used by solvePnP.
fit_corners_px = corners_px .+ 1.0
corner_results = Vector{Any}(undef, 4)
for corner_index in 1:4
    previous_corner_index = mod1(corner_index - 1, 4)
    next_corner_index = mod1(corner_index + 1, 4)
    corner = fit_corners_px[corner_index, :]
    previous_corner = fit_corners_px[previous_corner_index, :]
    next_corner = fit_corners_px[next_corner_index, :]

    edge_a = previous_corner .- corner
    edge_b = next_corner .- corner
    edge_a_length = norm(edge_a)
    edge_b_length = norm(edge_b)
    shortest_edge = min(edge_a_length, edge_b_length)
    max_profile_distance = edge_distance_fraction * shortest_edge
    max_profile_distance > profile_start_offset + min_profile_extent ||
        error("ArUco corner $corner_index is too close to its neighbours for the configured edge profiles.")

    corner_results[corner_index] = Camera.IntrinsicUncertainty.fit_corner(
        intensity_interpolant,
        corner,
        edge_a ./ edge_a_length,
        edge_b ./ edge_b_length;
        profile_positions=collect(range(
            profile_start_offset,
            max_profile_distance;
            length=n_profiles,
        )),
        half=min(40.0, shortest_edge / 12),
    )
    corner_results[corner_index] === nothing &&
        error("No usable sigmoid edge fit for ArUco corner $corner_index.")
end

Σ_corners = Matrix{Float64}[
    Matrix{Float64}(result.Σ) for result in corner_results
]
fitted_corners_px = copy(corners_px)
for corner_index in 1:4
    fitted_corners_px[corner_index, :] .+= corner_results[corner_index].q
end

results_path = joinpath(out_dir, "aruco_corner_results_$(zoom_level).jld2")
jldsave(
    results_path;
    corner_results,
    opencv_corners_px=corners_px,
    fitted_corners_px,
    Σ_corners,
    marker_ids,
    zoom_level,
)
println("Corner fits and covariances saved to: $results_path")

# ── plot corner corrections ───────────────────────────────────────────────────

correction_plot =
    Camera.IntrinsicUncertainty.plot_aruco_corner_corrections(corner_results)
correction_plot_path =
    joinpath(out_dir, "aruco_corner_corrections_$(zoom_level).png")
savefig(correction_plot, correction_plot_path)
display(correction_plot)

# ── save one sigmoid/profile diagnostic plot per corner ────────────────────────

for corner_index in eachindex(corner_results)
    result = corner_results[corner_index]
    result === nothing &&
        error("No usable corner fit for corner $corner_index.")

    diagnostic_plot =
        Camera.IntrinsicUncertainty.plot_aruco_corner_diagnostic(
            intensity,
            fit_corners_px[corner_index, :],
            result;
            corner_index=corner_index,
            edge_index=edge_sel,
            profile_index=profile_sel,
        )

    diagnostic_path = joinpath(
        out_dir,
        "aruco_corner_fit_$(zoom_level)_corner$(corner_index).png",
    )
    savefig(diagnostic_plot, diagnostic_path)
    display(diagnostic_plot)
end

# ── save one residual plot per corner ─────────────────────────────────────────

for corner_index in eachindex(corner_results)
    result = corner_results[corner_index]
    result === nothing &&
        error("No usable corner fit for corner $corner_index.")

    residual_plot =
        Camera.IntrinsicUncertainty.plot_aruco_edge_profile_residuals(
            intensity,
            result;
            corner_index=corner_index,
            edge_index=edge_sel,
        )

    residual_path = joinpath(
        out_dir,
        "aruco_edge_residuals_$(zoom_level)_corner$(corner_index)_edge$(edge_sel).png",
    )
    savefig(residual_plot, residual_path)
    display(residual_plot)
end

# ── bootstrap extrinsics ───────────────────────────────────────────────────────

# changed this to look at the cache from fit_corner 
boot_cache = joinpath(
    out_dir,
    "extrinsic_bootstrap_$(zoom_level)_fit_corner_N$(N_boot).jld2",
)
if isfile(boot_cache)
    cached_bootstrap = JLD2.load(boot_cache)
    required_cache_keys = (
        "R_samples", "t_samples", "N_boot", "Σ_corners", "seed",
        "zoom_level", "corners_px", "K", "dist", "marker_length_mm",
    )
    cache_inputs_match = all(haskey(cached_bootstrap, key) for key in required_cache_keys) &&
        cached_bootstrap["N_boot"] == N_boot &&
        cached_bootstrap["seed"] == seed &&
        cached_bootstrap["zoom_level"] == zoom_level &&
        cached_bootstrap["marker_length_mm"] == marker_length_mm &&
        isequal(cached_bootstrap["Σ_corners"], Σ_corners) &&
        isequal(cached_bootstrap["corners_px"], fitted_corners_px) &&
        isequal(cached_bootstrap["K"], K) &&
        isequal(cached_bootstrap["dist"], dist)
    if cache_inputs_match
        println("── Loading matching bootstrap cache: $boot_cache")
        R_samples = cached_bootstrap["R_samples"]
        t_samples = cached_bootstrap["t_samples"]
    else
        println("── Bootstrap cache inputs changed; recomputing...")
        R_samples, t_samples = Camera.IntrinsicUncertainty.bootstrap_extrinsics(
            fitted_corners_px, K, dist, marker_length_mm;
            Σ_corners=Σ_corners, N=N_boot, seed=seed,
        )
        jldsave(
            boot_cache;
            R_samples,
            t_samples,
            N_boot,
            Σ_corners,
            seed,
            zoom_level,
            corners_px=fitted_corners_px,
            K,
            dist,
            marker_length_mm,
        )
    end
else
    println("── No bootstrap cache found; starting bootstrap (N=$N_boot)...")
    R_samples, t_samples = Camera.IntrinsicUncertainty.bootstrap_extrinsics(
        fitted_corners_px, K, dist, marker_length_mm;
        Σ_corners=Σ_corners, N=N_boot, seed=seed,
    )
    jldsave(
        boot_cache;
        R_samples,
        t_samples,
        N_boot,
        Σ_corners,
        seed,
        zoom_level,
        corners_px=fitted_corners_px,
        K,
        dist,
        marker_length_mm,
    )
    println("Bootstrap saved: $boot_cache")
end

Σ_t, Σ_rvec, rvec_mean, t_mean =
    Camera.IntrinsicUncertainty.extrinsic_covariance(R_samples, t_samples)

# ── report uncertainty ────────────────────────────────────────────────────────

σ_rvec_deg = rad2deg.(sqrt.(diag(Σ_rvec)))
σ_t_mm = sqrt.(diag(Σ_t))
println("\nExtrinsic parameter standard deviations:")
for (name, value) in zip(("rx", "ry", "rz"), σ_rvec_deg)
    @printf("σ(%s) = %.6f deg\n", name, value)
end
for (name, value) in zip(("tx", "ty", "tz"), σ_t_mm)
    @printf("σ(%s) = %.6f mm\n", name, value)
end

# ── plot and save results ─────────────────────────────────────────────────────

plt = Camera.IntrinsicUncertainty.plot_extrinsic_distributions(
    R_samples, t_samples;
    save_fig=true, out_dir=out_dir, name_suffix=zoom_level,
)
display(plt)

plt_conv_mean = Camera.IntrinsicUncertainty.plot_extrinsic_convergence_mean(
    R_samples, t_samples;
    checkpoints = checkpoints_convergence,
    n_seeds     = n_seeds_convergence,
    save_fig    = true,
    out_dir     = out_dir,
    name_suffix = zoom_level,
)
display(plt_conv_mean)

plt_stab = Camera.IntrinsicUncertainty.plot_extrinsic_bootstrap_stability(
    R_samples, t_samples;
    checkpoints = checkpoints_stability,
    save_fig    = true,
    out_dir     = out_dir,
    name_suffix = zoom_level,
)
display(plt_stab)

result = (Σ_t, Σ_rvec, rvec_mean, t_mean)
output_path = joinpath(out_dir, "extrinsic_uncertainty_$(zoom_level).jld2")
Camera.IntrinsicUncertainty.save_extrinsic_uncertainty(result, output_path) ||
    error("Failed to save extrinsic uncertainty to $output_path")
println("\nExtrinsic uncertainty saved to: $output_path")
