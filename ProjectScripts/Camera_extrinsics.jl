project_root = normpath(joinpath(@__DIR__, ".."))
include(joinpath(project_root, "src", "camera.jl"))
using .Camera
using JLD2
using LinearAlgebra
using Printf

# ── configuration ─────────────────────────────────────────────────────────────

zoom_level = "zoom5"
intrinsics_file = joinpath(
    project_root, "data", "calibration_camera", zoom_level, "detected", "camint_all.jld2",
)
out_dir = joinpath(project_root, "data", "test_droplets")
img_file = joinpath(out_dir, "$(zoom_level).JPG")
marker_length_mm = 10.0

N_boot = 1000
checkpoints_convergence = unique(round.(Int,
    10 .^ range(log10(5), log10(N_boot), length=20)))
checkpoints_stability = filter(
    x -> x <= N_boot,
    [100, 200, 400, 800, 1500, 2500, 3500, 5000],
)
n_seeds_convergence = 10
Σ_corners = [
    [0.25 0.0; 0.0 0.25],  # corner 0: top-left [px²]
    [0.25 0.0; 0.0 0.25],  # corner 1: top-right [px²]
    [0.25 0.0; 0.0 0.25],  # corner 2: bottom-right [px²]
    [0.25 0.0; 0.0 0.25],  # corner 3: bottom-left [px²]
]
seed = 123

# ── load intrinsics ────────────────────────────────────────────────────────────

isfile(intrinsics_file) || error("Intrinsic calibration file not found: $intrinsics_file")
intrinsics = load(intrinsics_file)
K_raw = get(intrinsics, "camera_matrix", get(intrinsics, "K", nothing))
dist_raw = get(intrinsics, "dist_coeffs", get(intrinsics, "dist", nothing))
K_raw === nothing && error("No camera matrix found in $intrinsics_file")
dist_raw === nothing && error("No distortion coefficients found in $intrinsics_file")
K = Matrix{Float64}(K_raw)
dist = vec(Float64.(dist_raw))

# ── load or compute bootstrap samples ─────────────────────────────────────────

mkpath(out_dir)
boot_cache = joinpath(out_dir, "extrinsic_bootstrap_$(zoom_level)_N$(N_boot).jld2")
if isfile(boot_cache)
    println("── Bootstrap cache gevonden: $boot_cache")
    println("── Sla bootstrap over en laad bestaande samples...")
    cached_bootstrap = load(boot_cache)
    R_samples = cached_bootstrap["R_samples"]
    t_samples = cached_bootstrap["t_samples"]
else
    println("── Geen cache gevonden, start bootstrap (N=$N_boot)...")
    isfile(img_file) || error("Droplet image not found: $img_file")
    corners_px, marker_ids = Camera.IntrinsicUncertainty.detect_aruco_corners(
        img_file, K, dist,
    )
    println("Detected ArUco marker ID(s): ", join(marker_ids, ", "))

    R_samples, t_samples = Camera.IntrinsicUncertainty.bootstrap_extrinsics(
        corners_px, K, dist, marker_length_mm;
        Σ_corners=Σ_corners, N=N_boot, seed=seed,
    )
    jldsave(boot_cache; R_samples, t_samples, N_boot, Σ_corners, seed, zoom_level)
    println("── Bootstrap opgeslagen: $boot_cache")
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
