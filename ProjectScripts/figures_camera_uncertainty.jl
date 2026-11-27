using FileIO
using Images
using JLD2
using LinearAlgebra
using Plots
using Printf
using Random
using Statistics

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const RESULTS_FILE = joinpath(
    PROJECT_ROOT, "data", "test_droplets", "camtodrop_uncertainty",
    "uncertainty_camtodrop.jld2",
)
const INTRINSIC_BOOTSTRAP_FILE = joinpath(
    PROJECT_ROOT, "data", "calibration_camera", "zoom5", "detected",
    "subset_summaries", "reprojection_bootstrap_Nboot500.jld2",
)
const INTRINSIC_CALIBRATION_FILE = joinpath(
    PROJECT_ROOT, "data", "calibration_camera", "zoom5", "detected",
    "camint_all.jld2",
)
const ARUCO_IMAGE_FILE =
    joinpath(PROJECT_ROOT, "data", "test_droplets", "zoom5.JPG")
const FIGURE_DIR = joinpath(PROJECT_ROOT, "docs", "figures")
const INTRINSIC_LABELS = ["fx", "fy", "cx", "cy", "k1", "k2", "p1", "p2", "k3"]
const POSE_LABELS = ["δrx", "δry", "δrz", "δtx", "δty", "δtz"]
const PARAMETER_LABELS = vcat(INTRINSIC_LABELS, POSE_LABELS)

function require_keys(data, keys, path)
    missing = filter(key -> !haskey(data, key), keys)
    isempty(missing) ||
        error("$path is missing required keys: $(join(missing, ", ")).")
end

function correlation_from_covariance(covariance)
    standard_deviations = sqrt.(max.(diag(covariance), 0.0))
    denominator = standard_deviations * transpose(standard_deviations)
    return covariance ./ max.(denominator, eps(Float64))
end

function intrinsic_axis_values(samples, column, nominal)
    if column <= 4
        return 100 .* (samples[:, column] .- nominal[column]) ./ abs(nominal[column])
    end
    return samples[:, column] .- nominal[column]
end

function intrinsic_axis_label(column)
    return column <= 4 ? "Relative deviation (%)" :
           "Absolute coefficient deviation"
end

function parameter_axis_values(samples, column, nominal)
    if column <= 4
        return 100 .* (samples[:, column] .- nominal[column]) ./ abs(nominal[column])
    elseif column <= 9
        return samples[:, column] .- nominal[column]
    elseif column <= 12
        return 1_000 .* samples[:, column]
    end
    return 1_000 .* samples[:, column]
end

function parameter_axis_label(column)
    if column <= 4
        return "Relative deviation (%)"
    elseif column <= 9
        return "Absolute coefficient deviation"
    elseif column <= 12
        return "Rotation increment (mrad)"
    end
    return "Translation increment (μm)"
end

function save_figure(plot_object, filename)
    path = joinpath(FIGURE_DIR, filename)
    savefig(plot_object, path)
    println("Saved $path")
end

function strongest_intrinsic_pairs(correlation; count=4)
    pairs = [
        (abs(correlation[i, j]), i, j, correlation[i, j])
        for i in 1:8 for j in i+1:9
    ]
    sort!(pairs; by=first, rev=true)
    return pairs[1:min(count, length(pairs))]
end

function strongest_intrinsic_pose_pairs(correlation; count=6)
    pairs = [
        (abs(correlation[i, j]), i, j, correlation[i, j])
        for i in 1:9 for j in 10:15
    ]
    sort!(pairs; by=first, rev=true)
    return pairs[1:min(count, length(pairs))]
end

function covariance_ellipse(center, covariance; enlargement=1.0, points=160)
    decomposition = eigen(Symmetric(covariance))
    angles = range(0, 2π; length=points)
    unit_circle = vcat(reshape(cos.(angles), 1, :), reshape(sin.(angles), 1, :))
    axes = decomposition.vectors *
           Diagonal(sqrt.(max.(decomposition.values, 0.0))) * unit_circle
    return reshape(center, 2, 1) .+ enlargement .* axes
end

function contour_normals(r, z)
    n = length(r)
    normals = Matrix{Float64}(undef, n, 2)
    for i in 1:n
        left, right = i == 1 ? (1, 2) : i == n ? (n - 1, n) : (i - 1, i + 1)
        dr = r[right] - r[left]
        dz = z[right] - z[left]
        magnitude = hypot(dr, dz)
        magnitude > 0 || error("Degenerate contour tangent at point $i.")
        normals[i, :] .= (-dz / magnitude, dr / magnitude)
    end
    return normals
end

function normal_band_shape(r, z, normals, sigma; enlargement=10.0)
    upper_r = r .+ enlargement .* sigma .* normals[:, 1]
    upper_z = z .+ enlargement .* sigma .* normals[:, 2]
    lower_r = r .- enlargement .* sigma .* normals[:, 1]
    lower_z = z .- enlargement .* sigma .* normals[:, 2]
    return Shape(
        vcat(upper_r, reverse(lower_r)),
        vcat(upper_z, reverse(lower_z)),
    )
end

function mode_composition(covariance, r, z; count=5)
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
    scale = Vector{Float64}(undef, 2 * length(r))
    scale[1:2:end] .= r .- centroid[1]
    scale[2:2:end] .= z .- centroid[2]
    scale_residual = scale - Qrigid * (transpose(Qrigid) * scale)
    scale_residual ./= max(norm(scale_residual), eps(Float64))

    fractions = Matrix{Float64}(undef, min(count, length(order)), 3)
    for rank in axes(fractions, 1)
        vector = decomposition.vectors[:, order[rank]]
        rigid_fraction = sum(abs2, transpose(Qrigid) * vector)
        scale_fraction = abs2(dot(scale_residual, vector))
        fractions[rank, :] .= (
            rigid_fraction,
            scale_fraction,
            max(0.0, 1 - rigid_fraction - scale_fraction),
        )
    end
    return fractions
end

function main()
    for path in (
        RESULTS_FILE, INTRINSIC_BOOTSTRAP_FILE,
        INTRINSIC_CALIBRATION_FILE, ARUCO_IMAGE_FILE,
    )
        isfile(path) || error("Required saved input not found: $path")
    end
    mkpath(FIGURE_DIR)

    results = JLD2.load(RESULTS_FILE)
    bootstrap = JLD2.load(INTRINSIC_BOOTSTRAP_FILE)
    calibration = JLD2.load(INTRINSIC_CALIBRATION_FILE)
    require_keys(
        results,
        (
            "Σ_intrinsics", "Σ_ext_full", "Σ_camera_independent_full_ext",
            "Σ_camera_joint_reference", "Σ_corners", "fitted_aruco_corners_px",
            "virtual_r", "virtual_z", "normal_sigma_independent",
            "normal_sigma_joint", "normal_sigma_joint_camera", "normal_sigma_pixel",
            "normal_sigma_joint_total", "independent_apex_radius_mm",
            "joint_apex_radius_mm", "independent_parameter_samples",
            "joint_intrinsic_samples", "joint_pose_increments",
            "independent_contour_samples", "joint_contour_samples",
            "shared_corner_shift_assessment",
        ),
        RESULTS_FILE,
    )
    require_keys(
        bootstrap, ("samples", "mean_φ", "Σ_φ", "residuals", "Sigma_corner"),
        INTRINSIC_BOOTSTRAP_FILE,
    )

    Σ_intrinsics = Matrix{Float64}(results["Σ_intrinsics"])
    Σ_ext_full = Matrix{Float64}(results["Σ_ext_full"])
    φ_nominal = Float64[
        calibration["camera_matrix"][1, 1],
        calibration["camera_matrix"][2, 2],
        calibration["camera_matrix"][1, 3],
        calibration["camera_matrix"][2, 3],
        vec(calibration["dist_coeffs"])[1:5]...,
    ]
    φ_bootstrap_mean = Float64.(bootstrap["mean_φ"])
    intrinsic_samples = Matrix{Float64}(bootstrap["samples"])
    Σ_φ = Matrix{Float64}(bootstrap["Σ_φ"])
    ρ_φ = correlation_from_covariance(Σ_φ)

    p_relative = bar(
        INTRINSIC_LABELS[1:4],
        100 .* sqrt.(max.(diag(Σ_φ)[1:4], 0.0)) ./ abs.(φ_bootstrap_mean[1:4]);
        color=:steelblue, legend=false, ylabel="Relative σ (%)",
        title="Focal length and principal point",
    )
    p_distortion = bar(
        INTRINSIC_LABELS[5:9],
        sqrt.(max.(diag(Σ_φ)[5:9], 0.0));
        color=:darkorange, legend=false, ylabel="Absolute σ (coefficient)",
        title="Distortion coefficients",
    )
    save_figure(
        plot(p_relative, p_distortion; layout=(1, 2), size=(1050, 430)),
        "A1_intrinsic_marginal_sigma.pdf",
    )

    intrinsic_correlation_plot = heatmap(
        1:9, 1:9, ρ_φ;
        xticks=(1:9, INTRINSIC_LABELS),
        yticks=(1:9, INTRINSIC_LABELS),
        xrotation=45, yflip=true, aspect_ratio=:equal,
        clims=(-1, 1), c=cgrad([:royalblue, :white, :firebrick]),
        colorbar_title="Correlation",
        title="Intrinsic bootstrap correlation",
    )
    for i in 1:9, j in 1:9
        label_color = abs(ρ_φ[j, i]) > 0.62 ? :white : :black
        annotate!(
            intrinsic_correlation_plot, i, j,
            text(string(round(ρ_φ[j, i]; digits=2)), 7, label_color),
        )
    end
    save_figure(intrinsic_correlation_plot, "A2_intrinsic_correlation.pdf")

    intrinsic_pairs = strongest_intrinsic_pairs(ρ_φ; count=4)
    pair_plots = Plots.Plot[]
    for (_, i, j, correlation) in intrinsic_pairs
        x = intrinsic_axis_values(intrinsic_samples, i, φ_bootstrap_mean)
        y = intrinsic_axis_values(intrinsic_samples, j, φ_bootstrap_mean)
        push!(
            pair_plots,
            scatter(
                x, y;
                markersize=2.2, markerstrokewidth=0, alpha=0.5,
                xlabel="$(INTRINSIC_LABELS[i]) — $(intrinsic_axis_label(i))",
                ylabel="$(INTRINSIC_LABELS[j]) — $(intrinsic_axis_label(j))",
                title="$(INTRINSIC_LABELS[i])–$(INTRINSIC_LABELS[j]), ρ=$(round(correlation; digits=2))",
                label=false,
            ),
        )
    end
    save_figure(
        plot(pair_plots...; layout=(2, 2), size=(1000, 760)),
        "A3_intrinsic_pair_scatter.pdf",
    )

    image = FileIO.load(ARUCO_IMAGE_FILE)
    image_gray = Float64.(Gray.(image))
    corners = Matrix{Float64}(results["fitted_aruco_corners_px"]) .+ 1
    corner_covariances = results["Σ_corners"]
    enlargement = 20.0
    marker_span = max(
        maximum(corners[:, 1]) - minimum(corners[:, 1]),
        maximum(corners[:, 2]) - minimum(corners[:, 2]),
    )
    padding = max(40, round(Int, marker_span * 0.12))
    xlo = max(1, floor(Int, minimum(corners[:, 1]) - padding))
    xhi = min(size(image_gray, 2), ceil(Int, maximum(corners[:, 1]) + padding))
    ylo = max(1, floor(Int, minimum(corners[:, 2]) - padding))
    yhi = min(size(image_gray, 1), ceil(Int, maximum(corners[:, 2]) + padding))
    corner_plot = heatmap(
        xlo:xhi, ylo:yhi, image_gray[ylo:yhi, xlo:xhi];
        c=:grays, colorbar=false, yflip=true, aspect_ratio=:equal,
        xlabel="Image x (px)", ylabel="Image y (px)",
        title="Fitted ArUco corners and 2×2 covariance ellipses",
    )
    for i in 1:4
        next = mod1(i + 1, 4)
        plot!(
            corner_plot, corners[[i, next], 1], corners[[i, next], 2];
            color=:lime, linewidth=1.5, arrow=:closed,
            label=i == 1 ? "Marker edge directions" : false,
        )
    end
    for i in 1:4
        ellipse = covariance_ellipse(
            corners[i, :], Matrix{Float64}(corner_covariances[i]);
            enlargement,
        )
        plot!(
            corner_plot, ellipse[1, :], ellipse[2, :];
            color=:magenta, linewidth=3,
            label=i == 1 ? "Unit covariance ellipse ×$(Int(enlargement))" : false,
        )
        scatter!(
            corner_plot, [corners[i, 1]], [corners[i, 2]];
            color=:cyan, markershape=:xcross, markersize=7,
            label=i == 1 ? "Fitted corner" : false,
        )
        annotate!(corner_plot, corners[i, 1] + 8, corners[i, 2] + 8, text(string(i), 9, :yellow))
    end
    save_figure(corner_plot, "B1_aruco_corner_uncertainty.pdf")

    corner_noise_covariances = [
        Matrix{Float64}(covariance) for covariance in results["Σ_corners"]
    ]
    corner_noise_limit = 2.5 * maximum(
        sqrt(maximum(eigvals(Symmetric(covariance))))
        for covariance in corner_noise_covariances
    )
    corner_noise_plot = plot(
        layout=(2, 2), size=(1000, 820),
        xlabel="Δx (px)", ylabel="Δy (px)",
        xlims=(-corner_noise_limit, corner_noise_limit),
        ylims=(-corner_noise_limit, corner_noise_limit),
        aspect_ratio=:equal, legend=:bottomright,
    )
    corner_noise_rng = MersenneTwister(20261126)
    for i in 1:4
        covariance = corner_noise_covariances[i]
        decomposition = eigen(Symmetric(covariance))
        samples = decomposition.vectors *
                  Diagonal(sqrt.(max.(decomposition.values, 0.0))) *
                  randn(corner_noise_rng, 2, 500)
        plot_index = i
        scatter!(
            corner_noise_plot, samples[1, :], samples[2, :];
            subplot=plot_index, markersize=2, markerstrokewidth=0,
            alpha=0.35, color=:steelblue,
            label=i == 1 ? "Sample perturbations" : false,
        )
        for (sigma_level, style) in ((1.0, :solid), (2.0, :dash))
            ellipse = covariance_ellipse(zeros(2), covariance; enlargement=sigma_level)
            plot!(
                corner_noise_plot, ellipse[1, :], ellipse[2, :];
                subplot=plot_index, color=:black, linewidth=1.8,
                linestyle=style,
                label=i == 1 ? "$(Int(sigma_level))σ ellipse" : false,
            )
        end
        for (edge_index, neighbor) in enumerate((mod1(i - 1, 4), mod1(i + 1, 4)))
            edge_direction = corners[neighbor, :] .- corners[i, :]
            edge_direction ./= norm(edge_direction)
            endpoint = corner_noise_limit * 0.62 .* edge_direction
            plot!(
                corner_noise_plot, [0.0, endpoint[1]], [0.0, endpoint[2]];
                subplot=plot_index, color=edge_index == 1 ? :darkorange : :seagreen,
                linewidth=2, arrow=:closed,
                label=i == 1 ?
                      (edge_index == 1 ? "Edge to previous corner" : "Edge to next corner") :
                      false,
            )
        end
        annotate!(
            corner_noise_plot, -corner_noise_limit + 0.12,
            corner_noise_limit - 0.12,
            text(
                "σx=$(round(sqrt(covariance[1, 1]); digits=2)) px\n" *
                "σy=$(round(sqrt(covariance[2, 2]); digits=2)) px\n" *
                "ρ=$(round(correlation_from_covariance(covariance)[1, 2]; digits=3))",
                9, :left,
            );
            subplot=plot_index,
        )
        title!(corner_noise_plot, "Corner $i"; subplot=plot_index)
    end
    save_figure(corner_noise_plot, "B4_corner_noise_in_pixels.pdf")

    corner_covariance_8 = zeros(Float64, 8, 8)
    for i in 1:4
        indices = (2i - 1):(2i)
        corner_covariance_8[indices, indices] .= corner_noise_covariances[i]
    end
    corner_correlation_8 = correlation_from_covariance(corner_covariance_8)
    corner_coordinate_labels = ["x1", "y1", "x2", "y2", "x3", "y3", "x4", "y4"]
    corner_correlation_plot = heatmap(
        1:8, 1:8, corner_correlation_8;
        xticks=(1:8, corner_coordinate_labels),
        yticks=(1:8, corner_coordinate_labels),
        xrotation=45, yflip=true, aspect_ratio=:equal,
        clims=(-1, 1), c=cgrad([:royalblue, :white, :firebrick]),
        colorbar_title="Correlation",
        title="Corner-coordinate correlation (independent corners)",
        size=(800, 700),
    )
    for row in 1:8, column in 1:8
        annotate!(
            corner_correlation_plot, column, row,
            text(@sprintf("%.2f", corner_correlation_8[row, column]), 8, :black),
        )
    end
    for i in 1:4
        first_coordinate = 2i - 1
        last_coordinate = 2i
        plot!(
            corner_correlation_plot,
            [first_coordinate - 0.5, last_coordinate + 0.5,
             last_coordinate + 0.5, first_coordinate - 0.5,
             first_coordinate - 0.5],
            [first_coordinate - 0.5, first_coordinate - 0.5,
             last_coordinate + 0.5, last_coordinate + 0.5,
             first_coordinate - 0.5];
            color=:black, linewidth=1.5, label=false,
        )
    end
    save_figure(corner_correlation_plot, "B5_corner_correlation_matrix.pdf")

    shared_corner_assessment = results["shared_corner_shift_assessment"]
    shared_corner_labels = ["x shift", "y shift", "outward"]
    shared_corner_displacements = [
        case.median_abs_normal_displacement_um
        for case in shared_corner_assessment.cases
    ]
    shared_corner_plot = bar(
        shared_corner_labels, shared_corner_displacements;
        color=[:royalblue, :darkorange, :seagreen], legend=:topright,
        xlabel="Shared shift applied to all four corners",
        ylabel="Median |normal displacement| (µm)",
        title="Deterministic shared-corner sensitivity",
        label="Contour displacement",
    )
    hline!(
        shared_corner_plot,
        [shared_corner_assessment.median_joint_camera_sigma_um];
        color=:black, linestyle=:dash, linewidth=2,
        label="Median joint camera σ",
    )
    save_figure(shared_corner_plot, "E1_shared_corner_systematic_sensitivity.pdf")

    pose_sd = sqrt.(max.(diag(Σ_ext_full), 0.0))
    p_rotation = bar(
        POSE_LABELS[1:3], 1_000 .* pose_sd[1:3];
        color=:mediumpurple, legend=false, ylabel="σ (mrad)",
        title="Right-perturbation rotations",
    )
    p_translation = bar(
        POSE_LABELS[4:6], 1_000 .* pose_sd[4:6];
        color=:seagreen, legend=false, ylabel="σ (μm)",
        title="Pose translations",
    )
    save_figure(
        plot(p_rotation, p_translation; layout=(1, 2), size=(950, 420)),
        "B2_pose_marginal_sigma.pdf",
    )

    ρ_pose = correlation_from_covariance(Σ_ext_full)
    pose_correlation_plot = heatmap(
        1:6, 1:6, ρ_pose;
        xticks=(1:6, POSE_LABELS),
        yticks=(1:6, POSE_LABELS),
        xrotation=45, yflip=true, aspect_ratio=:equal,
        clims=(-1, 1), c=cgrad([:royalblue, :white, :firebrick]),
        colorbar_title="Correlation",
        title="Fixed-intrinsics pose bootstrap",
    )
    plot!(pose_correlation_plot, [3.5, 3.5], [0.5, 6.5]; color=:black, linewidth=1.4, label=false)
    plot!(pose_correlation_plot, [0.5, 6.5], [3.5, 3.5]; color=:black, linewidth=1.4, label=false)
    save_figure(pose_correlation_plot, "B3_pose_correlation.pdf")

    joint_parameter_samples = hcat(
        Matrix{Float64}(results["joint_intrinsic_samples"]),
        Matrix{Float64}(results["joint_pose_increments"]),
    )
    ρ_joint = correlation_from_covariance(cov(joint_parameter_samples; dims=1))
    Σ_independent_parameters = zeros(Float64, 15, 15)
    Σ_independent_parameters[1:9, 1:9] .= Σ_intrinsics
    Σ_independent_parameters[10:15, 10:15] .= Σ_ext_full
    ρ_independent = correlation_from_covariance(Σ_independent_parameters)

    function parameter_correlation_panel(correlation, title_text)
        panel = heatmap(
            1:15, 1:15, correlation;
            xticks=(1:15, PARAMETER_LABELS),
            yticks=(1:15, PARAMETER_LABELS),
            xrotation=60, yflip=true, aspect_ratio=:equal,
            clims=(-1, 1), c=cgrad([:royalblue, :white, :firebrick]),
            colorbar_title="ρ", title=title_text, tickfontsize=6,
        )
        plot!(panel, [9.5, 9.5], [0.5, 15.5]; color=:black, linewidth=1.5, label=false)
        plot!(panel, [0.5, 15.5], [9.5, 9.5]; color=:black, linewidth=1.5, label=false)
        return panel
    end
    p_joint = parameter_correlation_panel(ρ_joint, "Joint paired bootstrap")
    p_independent = parameter_correlation_panel(
        ρ_independent, "Independent model (cross block zero)",
    )
    save_figure(
        plot(p_joint, p_independent; layout=(1, 2), size=(1500, 700)),
        "C1_C2_joint_vs_independent_correlation.pdf",
    )

    coupled_pairs = strongest_intrinsic_pose_pairs(ρ_joint; count=6)
    independent_parameter_samples =
        Matrix{Float64}(results["independent_parameter_samples"])
    coupled_plots = Plots.Plot[]
    for (_, i, j, correlation) in coupled_pairs
        x_independent = parameter_axis_values(independent_parameter_samples, i, φ_nominal)
        y_independent = parameter_axis_values(independent_parameter_samples, j, φ_nominal)
        x_joint = parameter_axis_values(joint_parameter_samples, i, φ_nominal)
        y_joint = parameter_axis_values(joint_parameter_samples, j, φ_nominal)
        panel = scatter(
            x_independent, y_independent;
            markersize=1.8, markerstrokewidth=0, alpha=0.32,
            color=:gray35, label="Independent (500)",
            xlabel="$(PARAMETER_LABELS[i]) — $(parameter_axis_label(i))",
            ylabel="$(PARAMETER_LABELS[j]) — $(parameter_axis_label(j))",
            title="$(PARAMETER_LABELS[i])–$(PARAMETER_LABELS[j]), ρ=$(round(correlation; digits=2))",
        )
        scatter!(
            panel, x_joint, y_joint;
            markersize=1.6, markerstrokewidth=0, alpha=0.28,
            color=:darkorange, label="Joint (2000)",
        )
        push!(coupled_plots, panel)
    end
    save_figure(
        plot(coupled_plots...; layout=(3, 2), size=(1150, 1100)),
        "C3_intrinsic_pose_coupling.pdf",
    )

    r = Float64.(results["virtual_r"])
    z = Float64.(results["virtual_z"])
    contour_distance = vcat(0.0, cumsum(hypot.(diff(r), diff(z))))
    contour_fraction = 100 .* contour_distance ./ max(last(contour_distance), eps(Float64))
    p_sigma = plot(
        contour_fraction, 1_000 .* results["normal_sigma_independent"];
        label="Independent camera", linewidth=2,
        xlabel="Contour arclength (%)", ylabel="Normal σ (μm)",
        title="Normal-direction uncertainty",
    )
    plot!(p_sigma, contour_fraction, 1_000 .* results["normal_sigma_joint_camera"];
          label="Joint camera", linewidth=2)
    plot!(p_sigma, contour_fraction, 1_000 .* results["normal_sigma_pixel"];
          label="Pixel noise", linewidth=2)
    plot!(p_sigma, contour_fraction, 1_000 .* results["normal_sigma_joint_total"];
          label="Joint total", linewidth=2, linestyle=:dash)
    save_figure(p_sigma, "D1_normal_sigma_components.pdf")

    normals = contour_normals(r, z)
    band_enlargement = 10.0
    independent_band = normal_band_shape(
        1_000 .* r, 1_000 .* z, normals,
        1_000 .* Float64.(results["normal_sigma_independent"]);
        enlargement=band_enlargement,
    )
    joint_band = normal_band_shape(
        1_000 .* r, 1_000 .* z, normals,
        1_000 .* Float64.(results["normal_sigma_joint"]);
        enlargement=band_enlargement,
    )
    p_bands = plot(
        independent_band; color=:steelblue, fillalpha=0.18, linealpha=0,
        label="Independent 1σ ×$(Int(band_enlargement))",
        xlabel="r (μm)", ylabel="z (μm)", aspect_ratio=:equal,
        title="Virtual contour and enlarged normal bands",
    )
    plot!(p_bands, joint_band; color=:darkorange, fillalpha=0.24, linealpha=0,
          label="Joint 1σ ×$(Int(band_enlargement))")
    plot!(p_bands, 1_000 .* r, 1_000 .* z; color=:black, linewidth=1.5, label="Mean contour")
    save_figure(p_bands, "D2_contour_uncertainty_bands.pdf")

    independent_contours = Matrix{Float64}(results["independent_contour_samples"])
    joint_contours = Matrix{Float64}(results["joint_contour_samples"])
    sample_count = 30
    independent_indices =
        round.(Int, range(1, size(independent_contours, 1); length=sample_count))
    joint_indices =
        round.(Int, range(1, size(joint_contours, 1); length=sample_count))
    contour_enlargement = 10.0
    p_samples = plot(
        xlabel="r (μm)", ylabel="z (μm)", aspect_ratio=:equal,
        title="Sampled contours; deviations enlarged ×$(Int(contour_enlargement))",
    )
    for (sample_index, row) in enumerate(independent_indices)
        sample_r = view(independent_contours, row, 1:2:size(independent_contours, 2))
        sample_z = view(independent_contours, row, 2:2:size(independent_contours, 2))
        plotted_r = r .+ contour_enlargement .* (sample_r .- r)
        plotted_z = z .+ contour_enlargement .* (sample_z .- z)
        plot!(
            p_samples, 1_000 .* plotted_r, 1_000 .* plotted_z;
            color=:steelblue, linealpha=0.24, linewidth=0.7,
            label=sample_index == 1 ? "Independent (30 of 500)" : false,
        )
    end
    for (sample_index, row) in enumerate(joint_indices)
        sample_r = view(joint_contours, row, 1:2:size(joint_contours, 2))
        sample_z = view(joint_contours, row, 2:2:size(joint_contours, 2))
        plotted_r = r .+ contour_enlargement .* (sample_r .- r)
        plotted_z = z .+ contour_enlargement .* (sample_z .- z)
        plot!(
            p_samples, 1_000 .* plotted_r, 1_000 .* plotted_z;
            color=:darkorange, linealpha=0.24, linewidth=0.7,
            label=sample_index == 1 ? "Joint (30 of 2000)" : false,
        )
    end
    plot!(p_samples, 1_000 .* r, 1_000 .* z; color=:black, linewidth=1.6, label="Mean contour")
    save_figure(p_samples, "D3_sampled_contours.pdf")

    p_apex = histogram(
        1_000 .* results["independent_apex_radius_mm"];
        bins=32, normalize=:pdf, alpha=0.42, color=:steelblue,
        label="Independent (500)", xlabel="Apex-radius measure (μm)",
        ylabel="Probability density (1/μm)", title="Apex-radius measure",
    )
    histogram!(
        p_apex, 1_000 .* results["joint_apex_radius_mm"];
        bins=32, normalize=:pdf, alpha=0.42, color=:darkorange,
        label="Joint (2000)",
    )
    save_figure(p_apex, "D4_apex_radius_histogram.pdf")

    independent_modes = mode_composition(
        Matrix{Float64}(results["Σ_camera_independent_full_ext"]), r, z,
    )
    joint_modes = mode_composition(
        Matrix{Float64}(results["Σ_camera_joint_reference"]), r, z,
    )
    mode_labels = ["Mode $i" for i in 1:size(independent_modes, 1)]
    mode_positions = 1:length(mode_labels)
    mode_bar_width = 0.24
    p_independent_modes = bar(
        mode_positions .- mode_bar_width, 100 .* independent_modes[:, 1];
        bar_width=mode_bar_width, color=:steelblue, label="Rigid",
        xticks=(mode_positions, mode_labels), ylabel="Mode composition (%)",
        title="Independent camera",
    )
    bar!(
        p_independent_modes, mode_positions, 100 .* independent_modes[:, 2];
        bar_width=mode_bar_width, color=:darkorange, label="Scale",
    )
    bar!(
        p_independent_modes, mode_positions .+ mode_bar_width,
        100 .* independent_modes[:, 3];
        bar_width=mode_bar_width, color=:gray70, label="Shape",
    )
    p_joint_modes = bar(
        mode_positions .- mode_bar_width, 100 .* joint_modes[:, 1];
        bar_width=mode_bar_width, color=:steelblue, label="Rigid",
        xticks=(mode_positions, mode_labels), title="Joint camera",
    )
    bar!(
        p_joint_modes, mode_positions, 100 .* joint_modes[:, 2];
        bar_width=mode_bar_width, color=:darkorange, label="Scale",
    )
    bar!(
        p_joint_modes, mode_positions .+ mode_bar_width, 100 .* joint_modes[:, 3];
        bar_width=mode_bar_width, color=:gray70, label="Shape",
    )
    save_figure(
        plot(
            p_independent_modes, p_joint_modes;
            layout=(1, 2), size=(1250, 520), left_margin=14 * Plots.mm,
        ),
        "D5_covariance_modes.pdf",
    )
end

if abspath(PROGRAM_FILE) == (@__FILE__) || isinteractive()
    main()
end
