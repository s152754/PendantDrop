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
const BOOTSTRAP_FILE = joinpath(
    PROJECT_ROOT, "data", "calibration_camera", "zoom5", "detected",
    "subset_summaries", "reprojection_bootstrap_Nboot500.jld2",
)
const CALIBRATION_FILE = joinpath(
    PROJECT_ROOT, "data", "calibration_camera", "zoom5", "detected",
    "camint_all.jld2",
)
const MARKER_IMAGE = joinpath(PROJECT_ROOT, "data", "test_droplets", "zoom5.JPG")
const FIGURE_DIR = joinpath(PROJECT_ROOT, "docs", "figures", "student")
const MARKER_LENGTH_MM = 10.0
const CLOUD_SEED = 20261126

function save_student_figure(plot_object, filename)
    mkpath(FIGURE_DIR)
    path = joinpath(FIGURE_DIR, filename)
    savefig(plot_object, path)
    println("Saved $path")
end

function covariance_ellipse(covariance; sigma=1.0, points=160)
    decomposition = eigen(Symmetric(covariance))
    angles = range(0, 2π; length=points)
    unit_circle = vcat(reshape(cos.(angles), 1, :), reshape(sin.(angles), 1, :))
    return sigma .* decomposition.vectors *
           Diagonal(sqrt.(max.(decomposition.values, 0.0))) * unit_circle
end

function project_marker_corners(rotation, translation, K, dist)
    marker = Float64[
        -MARKER_LENGTH_MM / 2 MARKER_LENGTH_MM / 2 MARKER_LENGTH_MM / 2 -MARKER_LENGTH_MM / 2
         MARKER_LENGTH_MM / 2 MARKER_LENGTH_MM / 2 -MARKER_LENGTH_MM / 2 -MARKER_LENGTH_MM / 2
         0 0 0 0
    ]
    camera_points = rotation * marker .+ reshape(translation, 3, 1)
    x = camera_points[1, :] ./ camera_points[3, :]
    y = camera_points[2, :] ./ camera_points[3, :]
    r2 = x.^2 .+ y.^2
    radial = 1 .+ dist[1] .* r2 .+ dist[2] .* r2.^2 .+ dist[5] .* r2.^3
    x_distorted = x .* radial .+ 2dist[3] .* x .* y .+
                  dist[4] .* (r2 .+ 2x.^2)
    y_distorted = y .* radial .+ dist[3] .* (r2 .+ 2y.^2) .+
                  2dist[4] .* x .* y
    image_x = K[1, 1] .* x_distorted .+ K[1, 3] .+ 1
    image_y = K[2, 2] .* y_distorted .+ K[2, 3] .+ 1
    return hcat(image_x, image_y)
end

function amplified_contours(samples, reference_r, reference_z, indices, factor)
    return [
        (
            reference_r .+ factor .* (
                view(samples, row, 1:2:size(samples, 2)) .- reference_r
            ),
            reference_z .+ factor .* (
                view(samples, row, 2:2:size(samples, 2)) .- reference_z
            ),
        )
        for row in indices
    ]
end

function main()
    for path in (RESULTS_FILE, BOOTSTRAP_FILE, CALIBRATION_FILE, MARKER_IMAGE)
        isfile(path) || error("Required saved input not found: $path")
    end
    results = JLD2.load(RESULTS_FILE)
    bootstrap = JLD2.load(BOOTSTRAP_FILE)
    calibration = JLD2.load(CALIBRATION_FILE)
    K = Matrix{Float64}(calibration["camera_matrix"])
    dist = vec(Float64.(calibration["dist_coeffs"]))[1:5]
    intrinsics = Float64[
        K[1, 1], K[2, 2], K[1, 3], K[2, 3], dist...
    ]
    corner_covariances = [
        Matrix{Float64}(covariance) for covariance in results["Σ_corners"]
    ]
    corners = Matrix{Float64}(results["fitted_aruco_corners_px"])
    position_assessment = results["metadata"]["relative_position_assessment"]
    shared_assessment = results["shared_corner_shift_assessment"]

    pipeline = plot(
        xlims=(0, 12), ylims=(0, 3.2), aspect_ratio=:equal,
        axis=false, grid=false, legend=false, size=(1200, 360),
        title="A camera model turns pixels into a physical drop shape",
        titlefontsize=16,
    )
    stages = [
        ("Drop photo", "lightsteelblue"),
        ("Edge pixels", "lightcyan"),
        ("Camera model", "wheat"),
        ("Shape in mm", "honeydew"),
        ("Surface tension", "mistyrose"),
    ]
    centers = [1.2, 3.6, 6.0, 8.4, 10.8]
    for (index, (label, fill)) in enumerate(stages)
        x = centers[index]
        plot!(
            pipeline,
            Shape([x - 0.85, x + 0.85, x + 0.85, x - 0.85],
                  [1.05, 1.05, 2.05, 2.05]);
            fillcolor=fill, linecolor=:slategray, linewidth=1.5,
        )
        annotate!(pipeline, x, 1.55, text(label, 12, :black))
        if index < length(stages)
            plot!(
                pipeline, [x + 0.9, centers[index + 1] - 0.9], [1.55, 1.55];
                color=:slategray, linewidth=2, arrow=:closed,
            )
        end
    end
    annotate!(pipeline, 3.6, 0.55, text("find the outline", 10, :gray35))
    annotate!(pipeline, 6.0, 0.55, text("convert pixels using intrinsics + pose", 10, :gray35))
    annotate!(pipeline, 8.4, 0.55, text("measure the drop shape", 10, :gray35))
    save_student_figure(pipeline, "01_pixels_to_surface_tension.pdf")

    setup = plot(
        xlims=(-1.5, 6.4), ylims=(-2.2, 2.5), aspect_ratio=:equal,
        legend=false, grid=false, xlabel="Top view (schematic, not to scale)",
        yticks=false, size=(1000, 560),
        title="The marker tells us where the camera is relative to the drop",
    )
    plot!(
        setup, Shape([-0.65, 0.35, 0.35, -0.65], [-0.45, -0.45, 0.45, 0.45]);
        fillcolor=:lightsteelblue, linecolor=:slategray, linewidth=2,
    )
    plot!(setup, [0.35, 0.75], [0, 0]; color=:slategray, linewidth=5)
    annotate!(setup, -0.15, 0, text("Camera", 12, :black))
    plot!(setup, [0.75, 4.0], [0, 0.45]; color=:gray50, linestyle=:dash, linewidth=2)
    plot!(setup, [0.75, 4.0], [0, -0.85]; color=:gray50, linestyle=:dash, linewidth=2)
    plot!(
        setup, Shape([3.62, 4.38, 4.38, 3.62], [0.07, 0.07, 0.83, 0.83]);
        fillcolor=:white, linecolor=:black, linewidth=2,
    )
    for offset in (-0.25, 0.0, 0.25)
        plot!(setup, [4.0 + offset, 4.0 + offset], [0.08, 0.82];
              color=:gray50, linewidth=0.6)
        plot!(setup, [3.63, 4.37], [0.45 + offset, 0.45 + offset];
              color=:gray50, linewidth=0.6)
    end
    scatter!(setup, [4.0], [-0.85]; markersize=15, color=:darkorange, label=false)
    annotate!(setup, 4.45, -0.85, text("Drop", 12, :black))
    annotate!(setup, 4.0, 1.0, text("ArUco marker", 12, :black))
    plot!(setup, [0.1, 0.1], [1.65, 0.5]; color=:royalblue, arrow=:closed, linewidth=2)
    annotate!(setup, 0.25, 1.85, text("Intrinsic: inside the camera", 12, :royalblue, :left))
    annotate!(setup, 4.25, 1.95,
              text("Extrinsic: camera position and tilt\nrelative to the marker", 12, :darkgreen, :left))
    plot!(setup, [0.0, 4.0], [-1.6, -1.6]; color=:black, linewidth=1.5,
          arrow=(:closed, :closed))
    annotate!(setup, 2.0, -1.85, text("about 40 mm", 12, :black))
    save_student_figure(setup, "02_camera_marker_drop_setup.pdf")

    board = plot(
        xlims=(0, 7), ylims=(0, 7), aspect_ratio=:equal,
        legend=false, grid=false, axis=false, size=(560, 500),
        title="Calibration corners shift a little",
    )
    for row in 0:5, column in 0:5
        shade = iseven(row + column) ? :white : :gray20
        plot!(
            board,
            Shape([column, column + 1, column + 1, column],
                  [row, row, row + 1, row + 1]);
            fillcolor=shade, linecolor=:gray55, linewidth=0.7,
        )
    end
    true_points = [
        (column, row) for row in 1:5 for column in 1:5
        if (row + column) % 2 == 0
    ]
    selected = true_points[[1, 3, 4, 6, 7, 9, 10, 12, 13]]
    jitter = Float64[
        0.10 -0.08; -0.12 0.06; 0.06 0.11;
        -0.05 -0.12; 0.13 0.03; -0.10 0.09;
        0.08 -0.04; -0.06 0.12; 0.04 -0.09
    ]
    for (index, (x, y)) in enumerate(selected)
        plot!(board, [x, x + jitter[index, 1]], [y, y + jitter[index, 2]];
              color=:firebrick, alpha=0.65, linewidth=1)
    end
    scatter!(board, first.(selected), last.(selected);
             color=:black, markersize=4, label="Detected corner")
    scatter!(
        board,
        first.(selected) .+ jitter[:, 1],
        last.(selected) .+ jitter[:, 2];
        color=:firebrick, markersize=4, markerstrokewidth=0,
        label="Small shifted version",
    )
    focal_samples = Float64.(bootstrap["samples"][:, 1])
    focal_nominal = Float64(bootstrap["mean_φ"][1])
    focal_histogram = histogram(
        focal_samples; bins=18, color=:steelblue, linecolor=:white,
        xlabel="Focal length fx (pixels)", ylabel="Number of calibrations",
        legend=false, title="500 saved calibrations give a range",
    )
    vline!(focal_histogram, [focal_nominal]; color=:firebrick, linewidth=2)
    save_student_figure(
        plot(board, focal_histogram; layout=(1, 2), size=(1150, 500)),
        "03_small_corner_shifts_bootstrap.pdf",
    )

    image = Float64.(Gray.(FileIO.load(MARKER_IMAGE)))
    corners_plot = corners .+ 1
    marker_span = max(
        maximum(corners_plot[:, 1]) - minimum(corners_plot[:, 1]),
        maximum(corners_plot[:, 2]) - minimum(corners_plot[:, 2]),
    )
    padding = max(55, round(Int, marker_span * 0.10))
    xlo = max(1, floor(Int, minimum(corners_plot[:, 1]) - padding))
    xhi = min(size(image, 2), ceil(Int, maximum(corners_plot[:, 1]) + padding))
    ylo = max(1, floor(Int, minimum(corners_plot[:, 2]) - padding))
    yhi = min(size(image, 1), ceil(Int, maximum(corners_plot[:, 2]) + padding))
    marker_plot = heatmap(
        xlo:xhi, ylo:yhi, image[ylo:yhi, xlo:xhi];
        c=:grays, colorbar=false, yflip=true, aspect_ratio=:equal,
        xlabel="Image x (pixels)", ylabel="Image y (pixels)",
        title="Small corner errors move the fitted marker outline",
        size=(950, 700), legend=false,
    )
    pose_rotations = results["baseline_R_samples"]
    pose_translations = results["baseline_t_samples"]
    outline_indices = round.(Int, range(1, length(pose_rotations); length=30))
    for (draw_index, sample_index) in enumerate(outline_indices)
        projected = project_marker_corners(
            Matrix{Float64}(pose_rotations[sample_index]),
            Float64.(pose_translations[sample_index]),
            K, dist,
        )
        order = [1, 2, 3, 4, 1]
        plot!(
            marker_plot, projected[order, 1], projected[order, 2];
            color=:cyan, alpha=0.06, linewidth=0.7,
            label=draw_index == 1 ? "30 saved pose draws" : false,
        )
    end
    for i in 1:4
        ellipse = covariance_ellipse(corner_covariances[i]; sigma=20.0)
        plot!(
            marker_plot, corners_plot[i, 1] .+ ellipse[1, :],
            corners_plot[i, 2] .+ ellipse[2, :];
            color=:magenta, linewidth=2,
            label=i == 1 ? "Corner ellipse ×20" : false,
        )
        scatter!(
            marker_plot, [corners_plot[i, 1]], [corners_plot[i, 2]];
            color=:yellow, markershape=:xcross, markersize=7,
            label=i == 1 ? "Fitted corner" : false,
        )
    end
    save_student_figure(marker_plot, "04_marker_pose_bootstrap.pdf")

    first_covariance = corner_covariances[1]
    first_decomposition = eigen(Symmetric(first_covariance))
    corner_limit = 2.5 * sqrt(maximum(first_decomposition.values))
    corner_rng = MersenneTwister(CLOUD_SEED)
    corner_cloud = first_decomposition.vectors *
                   Diagonal(sqrt.(max.(first_decomposition.values, 0.0))) *
                   randn(corner_rng, 2, 500)
    corner_noise_plot = plot(
        xlims=(-corner_limit, corner_limit), ylims=(-corner_limit, corner_limit),
        aspect_ratio=:equal, xlabel="Δx (pixels)", ylabel="Δy (pixels)",
        title="One corner: x and y errors are nearly unrelated",
        legend=:bottomright, size=(650, 540),
    )
    scatter!(
        corner_noise_plot, corner_cloud[1, :], corner_cloud[2, :];
        color=:steelblue, alpha=0.35, markersize=2, markerstrokewidth=0,
        label="500 plotting draws",
    )
    ellipse = covariance_ellipse(first_covariance; sigma=1.0)
    plot!(
        corner_noise_plot, ellipse[1, :], ellipse[2, :];
        color=:black, linewidth=2, label="1σ ellipse",
    )
    first_edges = (
        corners_plot[4, :] .- corners_plot[1, :],
        corners_plot[2, :] .- corners_plot[1, :],
    )
    for (edge_index, direction) in enumerate(first_edges)
        direction ./= norm(direction)
        endpoint = 0.75 * corner_limit .* direction
        plot!(
            corner_noise_plot, [0.0, endpoint[1]], [0.0, endpoint[2]];
            color=edge_index == 1 ? :darkorange : :seagreen,
            linewidth=2, arrow=:closed,
            label=edge_index == 1 ? "Edge to previous corner" : "Edge to next corner",
        )
    end
    annotate!(
        corner_noise_plot, -corner_limit + 0.12, corner_limit - 0.9,
        text(
            "σx=$(round(sqrt(first_covariance[1, 1]); digits=2)) px\n" *
            "σy=$(round(sqrt(first_covariance[2, 2]); digits=2)) px\n" *
            "ρ=$(round(first_covariance[1, 2] /
                sqrt(first_covariance[1, 1] * first_covariance[2, 2]); digits=3))",
            10, :left,
        ),
    )

    edge_windows = plot(
        xlims=(-0.15, 1.15), ylims=(-0.15, 1.15), aspect_ratio=:equal,
        axis=false, grid=false, legend=false,
        title="Each corner uses its own stretch of each edge",
        size=(650, 540),
    )
    vertices = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    corner_colors = [:royalblue, :darkorange, :seagreen, :purple]
    boundary_order = [1, 2, 3, 4, 1]
    plot!(
        edge_windows, first.(vertices[boundary_order]),
        last.(vertices[boundary_order]);
        color=:black, linewidth=1.4,
    )
    for edge_index in 1:4
        next = mod1(edge_index + 1, 4)
        start_point, end_point = vertices[edge_index], vertices[next]
        plot!(
            edge_windows,
            [start_point[1] + 0.02 * (end_point[1] - start_point[1]),
             start_point[1] + 0.40 * (end_point[1] - start_point[1])],
            [start_point[2] + 0.02 * (end_point[2] - start_point[2]),
             start_point[2] + 0.40 * (end_point[2] - start_point[2])];
            color=corner_colors[edge_index], linewidth=7,
        )
        plot!(
            edge_windows,
            [end_point[1] + 0.02 * (start_point[1] - end_point[1]),
             end_point[1] + 0.40 * (start_point[1] - end_point[1])],
            [end_point[2] + 0.02 * (start_point[2] - end_point[2]),
             end_point[2] + 0.40 * (start_point[2] - end_point[2])];
            color=corner_colors[next], linewidth=7,
        )
    end
    annotate!(edge_windows, 0.5, 0.5, text("unused gap", 10, :gray35))
    corner_assumptions = plot(
        corner_noise_plot, edge_windows;
        layout=(1, 2), size=(1300, 560),
        plot_title="Separate edge windows support independent random fit errors",
    )
    save_student_figure(corner_assumptions, "05_corner_noise_assumptions.pdf")

    shared_schematic = plot(
        xlims=(-1.2, 4.2), ylims=(-1.1, 1.5), aspect_ratio=:equal,
        axis=false, grid=false, legend=false, size=(650, 500),
        title="A shared edge bias can make the marker look bigger",
    )
    plot!(
        shared_schematic, Shape([0, 1, 1, 0], [0, 0, 1, 1]);
        fillalpha=0, linecolor=:black, linewidth=2,
    )
    plot!(
        shared_schematic, Shape([-0.08, 1.08, 1.08, -0.08],
                                [-0.08, -0.08, 1.08, 1.08]);
        fillalpha=0, linecolor=:darkorange, linewidth=2, linestyle=:dash,
    )
    for (x, y) in ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
        dx = x == 0 ? -0.08 : 0.08
        dy = y == 0 ? -0.08 : 0.08
        plot!(
            shared_schematic, [x, x + dx], [y, y + dy];
            color=:darkorange, linewidth=1.5, arrow=:closed,
        )
    end
    annotate!(shared_schematic, 2.0, 0.72,
              text("same image size + bigger assumed marker\n→ camera seems closer", 11, :left))
    plot!(shared_schematic, [1.35, 3.2], [-0.65, -0.65];
          color=:slategray, linestyle=:dash, arrow=:closed, linewidth=1.5)
    annotate!(shared_schematic, 2.3, -0.9, text("estimated distance decreases", 10, :gray35))
    shared_values = [
        case.median_abs_normal_displacement_um for case in shared_assessment.cases
    ]
    shared_plot = bar(
        ["x +0.1 px", "y +0.1 px", "outward +0.1 px"],
        shared_values; color=[:royalblue, :darkorange, :seagreen],
        legend=false, ylabel="Median contour shift (µm)",
        title="These shared shifts are below camera uncertainty",
    )
    hline!(
        shared_plot, [shared_assessment.median_joint_camera_sigma_um];
        color=:black, linestyle=:dash, linewidth=2, label=false,
    )
    annotate!(
        shared_plot, 2.2, shared_assessment.median_joint_camera_sigma_um + 0.25,
        text("joint camera σ ≈ $(round(shared_assessment.median_joint_camera_sigma_um; digits=1)) µm", 9),
    )
    save_student_figure(
        plot(shared_schematic, shared_plot; layout=(1, 2), size=(1300, 560)),
        "06_shared_corner_errors.pdf",
    )

    independent_samples = Matrix{Float64}(results["independent_contour_samples"])
    joint_samples = Matrix{Float64}(results["joint_contour_samples"])
    reference_r = Float64.(results["virtual_r"])
    reference_z = Float64.(results["virtual_z"])
    independent_indices =
        round.(Int, range(1, size(independent_samples, 1); length=30))
    joint_indices =
        round.(Int, range(1, size(joint_samples, 1); length=30))
    exaggeration = 10.0
    independent_curves = amplified_contours(
        independent_samples, reference_r, reference_z,
        independent_indices, exaggeration,
    )
    joint_curves = amplified_contours(
        joint_samples, reference_r, reference_z, joint_indices, exaggeration,
    )
    all_r = vcat(
        reference_r,
        [curve[1] for curve in independent_curves]...,
        [curve[1] for curve in joint_curves]...,
    )
    all_z = vcat(
        reference_z,
        [curve[2] for curve in independent_curves]...,
        [curve[2] for curve in joint_curves]...,
    )
    r_padding = 0.04 * max(maximum(all_r) - minimum(all_r), eps(Float64))
    z_padding = 0.04 * max(maximum(all_z) - minimum(all_z), eps(Float64))
    contour_limits = (
        (minimum(all_r) - r_padding, maximum(all_r) + r_padding),
        (minimum(all_z) - z_padding, maximum(all_z) + z_padding),
    )
    function draw_contour_panel(curves, title_text, color)
        panel = plot(
            xlims=contour_limits[1], ylims=contour_limits[2],
            aspect_ratio=:equal, xlabel="Drop radius (mm)", ylabel="Drop height (mm)",
            title=title_text, legend=false,
        )
        for (index, (curve_r, curve_z)) in enumerate(curves)
            plot!(
                panel, curve_r, curve_z; color, alpha=0.28, linewidth=0.8,
                label=false,
            )
        end
        plot!(panel, reference_r, reference_z; color=:black, linewidth=1.5)
        return panel
    end
    independent_panel = draw_contour_panel(
        independent_curves, "Independent: wider spread", :steelblue,
    )
    joint_panel = draw_contour_panel(
        joint_curves, "Joint: errors partly cancel", :darkorange,
    )
    annotate!(independent_panel, contour_limits[1][1] + r_padding,
              contour_limits[2][2] - z_padding,
              text("30 outlines; deviations ×10", 9, :left))
    annotate!(joint_panel, contour_limits[1][1] + r_padding,
              contour_limits[2][2] - z_padding,
              text("30 outlines; deviations ×10", 9, :left))
    save_student_figure(
        plot(independent_panel, joint_panel; layout=(1, 2), size=(1250, 560)),
        "07_independent_vs_joint_contours.pdf",
    )

    independent_parameters =
        Matrix{Float64}(results["independent_parameter_samples"])
    joint_parameters = hcat(
        Matrix{Float64}(results["joint_intrinsic_samples"]),
        Matrix{Float64}(results["joint_pose_increments"]),
    )
    cy_nominal = intrinsics[4]
    cy_independent = independent_parameters[:, 4] .- cy_nominal
    ty_independent = 1_000 .* independent_parameters[:, 14]
    cy_joint = joint_parameters[:, 4] .- cy_nominal
    ty_joint = 1_000 .* joint_parameters[:, 14]
    coupling = cor(cy_joint, ty_joint)
    x_limits = extrema(vcat(cy_independent, cy_joint))
    y_limits = extrema(vcat(ty_independent, ty_joint))
    x_pad = 0.06 * max(x_limits[2] - x_limits[1], eps(Float64))
    y_pad = 0.06 * max(y_limits[2] - y_limits[1], eps(Float64))
    x_limits = (x_limits[1] - x_pad, x_limits[2] + x_pad)
    y_limits = (y_limits[1] - y_pad, y_limits[2] + y_pad)
    cy_ty_independent = scatter(
        cy_independent, ty_independent;
        xlims=x_limits, ylims=y_limits, color=:gray40,
        markersize=2, markerstrokewidth=0, alpha=0.35, legend=false,
        xlabel="cy (px)",
        ylabel="δty (µm)",
        guidefontsize=10, tickfontsize=9,
        bottom_margin=10Plots.mm,
        title="Independent draws: no cross-link",
    )
    cy_ty_joint = scatter(
        cy_joint, ty_joint;
        xlims=x_limits, ylims=y_limits, color=:darkorange,
        markersize=2, markerstrokewidth=0, alpha=0.35, legend=false,
        xlabel="cy (px)",
        ylabel="δty (µm)",
        guidefontsize=10, tickfontsize=9,
        bottom_margin=10Plots.mm,
        title="Joint draws: opposite changes pair up",
    )
    annotate!(
        cy_ty_joint, x_limits[1] + x_pad, y_limits[2] - y_pad,
        text("ρ=$(round(coupling; digits=2))", 11, :left),
    )
    save_student_figure(
        plot(
            cy_ty_independent, cy_ty_joint;
            layout=(1, 2), size=(1250, 600),
            left_margin=14 * Plots.mm, bottom_margin=10Plots.mm,
        ),
        "08_strongest_joint_pair.pdf",
    )

    offset_sigma_um = 1_000 * position_assessment.break_even_offset_mm *
                     position_assessment.median_normal_sigma_per_mm_offset
    summary_values = 1_000 .* [
        median(results["normal_sigma_independent"]),
        median(results["normal_sigma_joint_camera"]),
        median(results["normal_sigma_pixel"]),
        offset_sigma_um / 1_000,
    ]
    offset_label = "Depth offset\n$(round(position_assessment.break_even_offset_mm; digits=2)) mm"
    summary_labels = ["Independent\ncamera", "Joint\ncamera", "Pixel noise\nassumed", offset_label]
    summary_plot = bar(
        summary_labels, summary_values;
        color=[:steelblue, :darkorange, :gray65, :seagreen],
        legend=false, ylabel="Median normal uncertainty (µm)",
        title="Depth offset matches joint camera uncertainty",
        titlefontsize=14,
        ylims=(0, 36),
        xtickfontsize=10, ytickfontsize=11,
    )
    for (index, value) in enumerate(summary_values)
        annotate!(summary_plot, index, value + 0.8, text("$(round(value; digits=1))", 10))
    end
    save_student_figure(summary_plot, "09_summary_uncertainty.pdf")
end

if abspath(PROGRAM_FILE) == (@__FILE__) || isinteractive()
    main()
end
