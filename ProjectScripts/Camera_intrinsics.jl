##  intrinsic_uncertainty.jl
##
##  Uncertainty in camera intrinsic parameters.
##
##  Sources (one section per source):
##    [1] Number of calibration images  ← implemented here
##    [2] Reprojection / corner noise   ← implemented
##    [3] Corner-detection bias         ← TODO

project_root = normpath(joinpath(@__DIR__, ".."))
include(joinpath(project_root, "src", "camera.jl"))
using .Camera
using Printf
using Statistics
using LinearAlgebra

# ── configuration ─────────────────────────────────────────────────────────────

zoom        = "zoom5"
img_dir     = joinpath(project_root, "data", "calibration_camera", zoom)
out_dir     = joinpath(img_dir, "detected")
summary_dir = joinpath(out_dir, "subset_summaries")

subset_sizes = [5]
n_seeds      = 1

# ── [1] uncertainty from number of calibration images ────────────────────────

# Matching subset_<size>_seeds<count>.jld2 files are loaded and reused for plotting.
println("\n── [1] Subset studie ──────────────────────────────")
println("Cache wordt automatisch hergebruikt indien aanwezig.")
t1 = time()
results, K_ref, dist_ref = Camera.IntrinsicUncertainty.run_subset_study(
    img_dir,
    out_dir,
    summary_dir;
    subset_sizes = subset_sizes,
    n_seeds      = n_seeds,
)
println("── [1] klaar in $(round(time()-t1, digits=1)) s ──")

summary = Camera.IntrinsicUncertainty.summarise_by_subset(results, subset_sizes)

p_intrinsics, p_distortion = Camera.IntrinsicUncertainty.plot_subset_influence(
    summary,
    K_ref,
    dist_ref;
    save_dir = summary_dir,
)

p_intrinsics_by_seed, p_distortion_by_seed = Camera.IntrinsicUncertainty.plot_subset_seed_values(
    results,
    K_ref,
    dist_ref;
    save_dir = summary_dir,
)

if p_intrinsics !== nothing
    display(p_intrinsics)
    display(p_distortion)
end

if p_intrinsics_by_seed !== nothing
    display(p_intrinsics_by_seed)
    display(p_distortion_by_seed)
end

# ── [2] uncertainty from reprojection error / corner noise ───────────────────

N_boot = 500
println("\n── [2] Bootstrap ──────────────────────────────────")
println("N_boot = $N_boot")
println("Cache wordt automatisch hergebruikt indien aanwezig.")
t2 = time()
sigma_x, sigma_y, residuals, Sigma_corner, samples, mean_φ, Σ_φ, N_boot =
    Camera.IntrinsicUncertainty.run_reprojection_bootstrap(
        out_dir,
        summary_dir;
        N_boot = N_boot,
    )
println("── [2] klaar in $(round(time()-t2, digits=1)) s ──")

println("Sigma_corner (px²):")
display(Sigma_corner)

p_residuals = Camera.IntrinsicUncertainty.plot_reprojection_residuals(
    residuals;
    save_dir = summary_dir,
)
display(p_residuals)

println("Covariance matrix Σ_φ:")
display(Σ_φ)

println("\n── Consistentiecheck: bootstrap σ vs std_OpenCV ──────────")
Printf.@printf("%-7s | %-12s | %-12s | %-8s\n",
               "param", "σ_bootstrap", "std_OpenCV", "ratio")
println(repeat("-", 50))

param_names_check = ["fx", "fy", "cx", "cy", "k1", "k2", "p1", "p2", "k3"]
ref_result = Camera.IntrinsicUncertainty.calibrate_camera_from_images(
    Camera.IntrinsicUncertainty.list_camera_images(img_dir),
)
get(ref_result, "success", false) ||
    error("Multi-image reference calibration failed for OpenCV standard deviations.")
std_opencv = ref_result["std_intrinsics"][1:9]
sigma_bootstrap = sqrt.(LinearAlgebra.diag(Σ_φ))

for (i, name) in enumerate(param_names_check)
    ratio = sigma_bootstrap[i] / std_opencv[i]
    flag = 0.5 < ratio < 2.0 ? "✓" : "⚠"
    Printf.@printf("%-7s | %12.6f | %12.6f | %6.2f %s\n",
                   name, sigma_bootstrap[i], std_opencv[i], ratio, flag)
end

println()
println("ratio ≈ 1.0 → bootstrap en OpenCV consistent")
println("ratio >> 1  → bootstrap schat meer onzekerheid")
println("ratio << 1  → bootstrap schat minder onzekerheid")

# ── [3] convergentieanalyse: gecombineerde OpenCV kalibratie per subset ───

subset_sizes_conv = [5, 10, 20, 35, 50]
seed_conv = 1
all_files = Camera.IntrinsicUncertainty.list_camera_images(img_dir)

println("\n── [3] Convergentieanalyse ────────────────────────────")
println("Aantal beschikbare beelden: $(length(all_files))")
println("Subset-groottes: $subset_sizes_conv")
println("Seed: $seed_conv\n")
conv_cache = joinpath(summary_dir, "convergence_subset.jld2")
if isfile(conv_cache)
    println("── [3] Laden uit cache: $conv_cache")
    loaded = Camera.IntrinsicUncertainty.load(conv_cache)
    conv_results = loaded["conv_results"]
else
    println("── [3] Berekenen...")
    t3 = time()
    conv_results = Dict{Int,Dict{String,Any}}()
    for N in subset_sizes_conv
        tN = time()
        println("Kalibratie met N=$N beelden...")
        selected_files = Camera.IntrinsicUncertainty.random_image_subset(
            img_dir;
            subset_size=N,
            seed=seed_conv,
        )
        result = Camera.IntrinsicUncertainty.calibrate_camera_from_images(selected_files)
        get(result, "success", false) ||
            error("Camera calibration failed for subset size N=$N")

        K = Camera.IntrinsicUncertainty._coerce_camera_matrix(result["camera_matrix"])
        dist = Camera.IntrinsicUncertainty._coerce_distortion_coeffs(result["dist_coeffs"])
        std_int = result["std_intrinsics"][1:9]
        conv_results[N] = Dict(
            "K" => K,
            "dist" => dist,
            "std_intrinsics" => std_int,
        )
        println("  ✓ klaar in $(round(time()-tN, digits=1)) s  " *
                "| fx=$(round(K[1,1], digits=1))  " *
                "| RMS=$(round(result["rms"], digits=3)) px")
    end
    Camera.IntrinsicUncertainty.jldsave(conv_cache; conv_results)
    println("── [3] klaar in $(round(time()-t3, digits=1)) s ──")
end

println("── Resultaten ─────────────────────────────────────────\n")

param_names = ["fx", "fy", "cx", "cy", "k1", "k2", "p1", "p2", "k3"]
K_ref_matrix = Camera.IntrinsicUncertainty._coerce_camera_matrix(K_ref)
dist_ref_vector = Camera.IntrinsicUncertainty._coerce_distortion_coeffs(dist_ref)
multi_values = vcat(
    Float64[K_ref_matrix[1, 1], K_ref_matrix[2, 2], K_ref_matrix[1, 3], K_ref_matrix[2, 3]],
    dist_ref_vector[1:5],
)

for (i, param_name) in enumerate(param_names)
    println("=== $param_name ===")
    Printf.@printf("%-7s | %-12s | %-12s\n", "N", "waarde", "std_OpenCV")
    for N in subset_sizes_conv
        entry = conv_results[N]
        K = entry["K"]
        dist = entry["dist"]
        values = vcat(
            Float64[K[1, 1], K[2, 2], K[1, 3], K[2, 3]],
            dist[1:5],
        )
        Printf.@printf("%-7d | %12.6f | %12.6f\n",
                       N, values[i], entry["std_intrinsics"][i])
    end
    Printf.@printf("%-7s | %12.6f | %12s\n", "multi", multi_values[i], "-")

    entry_50 = conv_results[50]
    K_50 = entry_50["K"]
    dist_50 = entry_50["dist"]
    values_50 = vcat(
        Float64[K_50[1, 1], K_50[2, 2], K_50[1, 3], K_50[2, 3]],
        dist_50[1:5],
    )
    if abs(values_50[i] - multi_values[i]) < 2 * entry_50["std_intrinsics"][i]
        println("✓ consistent met multi-image referentie")
    else
        println("⚠ afwijking van multi-image referentie")
    end
end

# ── [4] figuren ────────────────────────────────────────────────────────────

fig_dir = joinpath(summary_dir, "figures")
mkpath(fig_dir)
println("\n── [4] Figuren ────────────────────────────────────")

boot_cache = Camera.IntrinsicUncertainty.load(joinpath(
    summary_dir, "reprojection_bootstrap_Nboot$(N_boot).jld2"))
Sigma_corner_fig = Matrix{Float64}(boot_cache["Sigma_corner"])
samples_fig = Matrix{Float64}(boot_cache["samples"])
mean_φ_fig = Vector{Float64}(boot_cache["mean_φ"])
Σ_φ_fig = Matrix{Float64}(boot_cache["Σ_φ"])

if haskey(boot_cache, "residuals")
    residuals_raw = boot_cache["residuals"]
else
    println("  Residuals herberekenen...")
    ref = Camera.IntrinsicUncertainty.load(joinpath(out_dir, "camint_all.jld2"))
    K_r = Camera.IntrinsicUncertainty._coerce_camera_matrix(ref["camera_matrix"])
    d_r = Camera.IntrinsicUncertainty._coerce_distortion_coeffs(ref["dist_coeffs"])
    _, _, residuals_raw, _ = Camera.IntrinsicUncertainty.compute_reprojection_residuals(
        ref["objpoints"], ref["imgpoints"], ref["rvecs"], ref["tvecs"], K_r, d_r)
end
residuals_fig = residuals_raw isa AbstractMatrix ?
    Matrix{Float64}(residuals_raw) :
    reduce(vcat, permutedims.(residuals_raw))

P = Camera.IntrinsicUncertainty.Plots

# ── Figuur 1: reprojection residuals met Sigma_corner ellips ──
println("  Figuur 1: reprojection residuals...")
p1_fig = P.scatter(
    residuals_fig[:, 1], residuals_fig[:, 2];
    marker=:circle, markersize=2, alpha=0.4,
    xlabel="x residual [px]", ylabel="y residual [px]",
    title="Reprojection residuals with Sigma_corner ellipse",
    legend=:topright, label="residuals",
    aspect_ratio=:equal)

F = LinearAlgebra.eigen(Sigma_corner_fig)
θ = range(0, 2π, length=200)
for scale in [1.0, 2.0]
    ex = scale .* (sqrt(F.values[1]) .* cos.(θ) .* F.vectors[1, 1] .+
                   sqrt(F.values[2]) .* sin.(θ) .* F.vectors[1, 2])
    ey = scale .* (sqrt(F.values[1]) .* cos.(θ) .* F.vectors[2, 1] .+
                   sqrt(F.values[2]) .* sin.(θ) .* F.vectors[2, 2])
    P.plot!(p1_fig, ex, ey;
            linewidth=2, label="$(Int(scale))σ ellipse")
end
P.savefig(p1_fig, joinpath(fig_dir, "reprojection_residuals_ellipse.png"))
println("  ✓ saved: reprojection_residuals_ellipse.png")

# ── Figuur 2: convergentie intrinsieke parameters ──
println("  Figuur 2: convergentie intrinsics...")
conv_loaded = Camera.IntrinsicUncertainty.load(conv_cache)
conv_results = conv_loaded["conv_results"]
Ns = sort(collect(keys(conv_results)))
conv_params = Dict(
    N => vcat(
        Float64[
            conv_results[N]["K"][1, 1], conv_results[N]["K"][2, 2],
            conv_results[N]["K"][1, 3], conv_results[N]["K"][2, 3],
        ],
        conv_results[N]["dist"][1:5],
    )
    for N in Ns
)

intrinsic_names = ["fx", "fy", "cx", "cy"]
intrinsic_idx = [1, 2, 3, 4]
panels_int = []
for (name, idx) in zip(intrinsic_names, intrinsic_idx)
    vals = [conv_params[N][idx] for N in Ns]
    stds = [conv_results[N]["std_intrinsics"][idx] for N in Ns]
    p = P.plot(Ns, vals;
               ribbon=stds, fillalpha=0.3,
               marker=:circle, xlabel="N images",
               ylabel="value", title=name, legend=false)
    P.hline!(p, [vals[end]]; linestyle=:dash, color=:black)
    push!(panels_int, p)
end
p2_fig = P.plot(panels_int...;
                layout=(2, 2), size=(900, 700),
                plot_title="Convergence of intrinsic parameters")
P.savefig(p2_fig, joinpath(fig_dir, "convergence_intrinsics.png"))
println("  ✓ saved: convergence_intrinsics.png")

# ── Figuur 3: convergentie distorsie parameters ──
println("  Figuur 3: convergentie distorsie...")
distortion_names = ["k1", "k2", "p1", "p2", "k3"]
distortion_idx = [5, 6, 7, 8, 9]
panels_dist = []
for (name, idx) in zip(distortion_names, distortion_idx)
    vals = [conv_params[N][idx] for N in Ns]
    stds = [conv_results[N]["std_intrinsics"][idx] for N in Ns]
    p = P.plot(Ns, vals;
               ribbon=stds, fillalpha=0.3,
               marker=:circle, xlabel="N images",
               ylabel="value", title=name, legend=false)
    P.hline!(p, [vals[end]]; linestyle=:dash, color=:black)
    push!(panels_dist, p)
end
p3_fig = P.plot(panels_dist...;
                layout=(2, 3), size=(1100, 700),
                plot_title="Convergence of distortion parameters")
P.savefig(p3_fig, joinpath(fig_dir, "convergence_distortion.png"))
println("  ✓ saved: convergence_distortion.png")

# ── Figuur 4: bootstrap histogrammen ──
println("  Figuur 4: bootstrap histogrammen...")
param_names_boot = ["fx", "fy", "cx", "cy", "k1", "k2", "p1", "p2", "k3"]
panels_boot = []
for (i, name) in enumerate(param_names_boot)
    p = P.histogram(samples_fig[:, i];
                    bins=30, title=name, legend=false,
                    xlabel="value", ylabel="count")
    P.vline!(p, [mean_φ_fig[i]];
             linewidth=2, color=:red, label="mean")
    push!(panels_boot, p)
end
p4_fig = P.plot(panels_boot...;
                layout=(3, 3), size=(1200, 900),
                plot_title="Bootstrap parameter distributions (N_boot=$N_boot)")
P.savefig(p4_fig, joinpath(fig_dir, "bootstrap_histograms.png"))
println("  ✓ saved: bootstrap_histograms.png")

# ── Figuur 5: bootstrap convergentie ──
println("  Figuur 5: bootstrap convergentie...")

checkpoints = filter(x -> x <= N_boot, [10, 25, 50, 100, 150, 200, 300, 400, 500])
param_names_boot = ["fx", "fy", "cx", "cy", "k1", "k2", "p1", "p2", "k3"]

std_per_checkpoint = zeros(length(checkpoints), 9)
for (j, M) in enumerate(checkpoints)
    for i in 1:9
        std_per_checkpoint[j, i] = Statistics.std(samples_fig[1:M, i])
    end
end

std_final = [Statistics.std(samples_fig[:, i]) for i in 1:9]

panels_conv_int = []
for (i, name) in enumerate(["fx", "fy", "cx", "cy"])
    p = P.plot(checkpoints, std_per_checkpoint[:, i];
               marker=:circle, xlabel="N_boot samples",
               ylabel="std", title=name, legend=false)
    P.hline!(p, [std_final[i]];
             linestyle=:dash, color=:black)
    push!(panels_conv_int, p)
end
p5a = P.plot(panels_conv_int...;
             layout=(2, 2), size=(900, 700),
             plot_title="Bootstrap convergence — intrinsics")
P.savefig(p5a, joinpath(fig_dir, "bootstrap_convergence_intrinsics.png"))
println("  ✓ saved: bootstrap_convergence_intrinsics.png")

panels_conv_dist = []
for (i, name) in enumerate(["k1", "k2", "p1", "p2", "k3"])
    idx = i + 4
    p = P.plot(checkpoints, std_per_checkpoint[:, idx];
               marker=:circle, xlabel="N_boot samples",
               ylabel="std", title=name, legend=false)
    P.hline!(p, [std_final[idx]];
             linestyle=:dash, color=:black)
    push!(panels_conv_dist, p)
end
p5b = P.plot(panels_conv_dist...;
             layout=(2, 3), size=(1100, 700),
             plot_title="Bootstrap convergence — distortion")
P.savefig(p5b, joinpath(fig_dir, "bootstrap_convergence_distortion.png"))
println("  ✓ saved: bootstrap_convergence_distortion.png")

println("\n── [4] Alle figuren opgeslagen in: $fig_dir ──")
