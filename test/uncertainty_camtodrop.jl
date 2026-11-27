using ForwardDiff
using FileIO
using JLD2
using LinearAlgebra
using Test

include("../src/forwardmodel.jl")
using .Model

const UQ = Model.Camera.IntrinsicUncertainty
const CAM = Model.Camera

@testset "Camera-to-droplet uncertainty helpers" begin
    intrinsics = Float64[
        800.0, 810.0, 320.0, 240.0,
        0.01, -0.001, 0.0002, -0.0001, 0.0001,
    ]
    R_mean = UQ._right_perturbed_rotation(
        Matrix{Float64}(I, 3, 3), Float64[0.01, -0.02, 0.03],
    )
    t_mean = Float64[10.0, -3.0, 500.0]
    R_rel = CAM.getRrel_fn(3.0, -4.0, 7.0)
    t_rel = Float64[1.0, -0.5, 2.0]
    parameters = vcat(intrinsics, zeros(6), 0.0)
    r = Float64[-10.0, 0.0, 10.0]
    z = Float64[-5.0, 3.0, 8.0]

    u, v, _ = UQ.forward_project_camtodrop(
        r, z, parameters, R_mean, t_mean, R_rel, t_rel;
        plane_coordinate=0.0,
    )
    recovered_r, recovered_z = UQ.inverse_project_camtodrop(
        u, v, parameters, R_mean, t_mean, R_rel, t_rel;
        plane_coordinate=0.0,
    )
    @test recovered_r ≈ r atol=1e-8
    @test recovered_z ≈ z atol=1e-8

    forward_fx_derivative = ForwardDiff.derivative(intrinsics[1]) do fx
        varied = vcat(fx, parameters[2:end])
        first(UQ.forward_project_camtodrop(
            r, z, varied, R_mean, t_mean, R_rel, t_rel;
            plane_coordinate=0.0,
        ))[1]
    end
    @test isfinite(forward_fx_derivative)
    @test !iszero(forward_fx_derivative)

    jacobians = UQ.camtodrop_jacobians(
        u, v, intrinsics, R_mean, t_mean, R_rel, t_rel;
        plane_coordinate=0.0,
    )
    @test size(jacobians.J_intrinsics) == (2 * length(r), 9)
    @test size(jacobians.J_extrinsics) == (2 * length(r), 6)
    @test size(jacobians.J_offset) == (2 * length(r), 1)
    @test size(jacobians.J_pixels) == (2 * length(r), 2 * length(r))
    @test all(isfinite, jacobians.J_intrinsics)
    @test all(isfinite, jacobians.J_extrinsics)

    step = 1e-3
    varied_plus, varied_minus = copy(parameters), copy(parameters)
    varied_plus[1] += step
    varied_minus[1] -= step
    plus_r, plus_z = UQ.inverse_project_camtodrop(
        u, v, varied_plus, R_mean, t_mean, R_rel, t_rel;
    )
    minus_r, minus_z = UQ.inverse_project_camtodrop(
        u, v, varied_minus, R_mean, t_mean, R_rel, t_rel;
    )
    finite_difference_fx = (vcat(plus_r[1], plus_z[1]) .-
                            vcat(minus_r[1], minus_z[1])) ./ (2 * step)
    @test jacobians.J_intrinsics[[1, 2], 1] ≈ finite_difference_fx atol=1e-5

    for parameter_index in (2, 4)
        plus, minus = copy(parameters), copy(parameters)
        plus[parameter_index] += step
        minus[parameter_index] -= step
        plus_r, plus_z = UQ.inverse_project_camtodrop(
            u, v, plus, R_mean, t_mean, R_rel, t_rel,
        )
        minus_r, minus_z = UQ.inverse_project_camtodrop(
            u, v, minus, R_mean, t_mean, R_rel, t_rel,
        )
        finite_difference = (
            vcat(plus_r[1], plus_z[1]) .-
            vcat(minus_r[1], minus_z[1])
        ) ./ (2 * step)
        @test jacobians.J_intrinsics[[1, 2], parameter_index] ≈
              finite_difference atol=1e-5
    end

    pixel_plus = vcat(u[1], v[1])
    pixel_minus = copy(pixel_plus)
    pixel_plus[2] += 1e-3
    pixel_minus[2] -= 1e-3
    plus_r, plus_z = UQ.inverse_project_camtodrop(
        pixel_plus[1:1], pixel_plus[2:2], parameters, R_mean, t_mean, R_rel, t_rel,
    )
    minus_r, minus_z = UQ.inverse_project_camtodrop(
        pixel_minus[1:1], pixel_minus[2:2], parameters, R_mean, t_mean, R_rel, t_rel,
    )
    finite_difference_pixel = (
        vcat(plus_r[1], plus_z[1]) .-
        vcat(minus_r[1], minus_z[1])
    ) ./ 2e-3
    @test jacobians.J_pixels[[1, 2], 2] ≈ finite_difference_pixel atol=1e-5

    Σ_intrinsics = Matrix(Diagonal(fill(1e-6, 9)))
    Σ_ext = Matrix(Diagonal(Float64[1e-6, 1e-6, 1e-6, 1e-4, 1e-4, 1e-4]))
    covariance =
        jacobians.J_intrinsics * Σ_intrinsics * transpose(jacobians.J_intrinsics) +
        jacobians.J_extrinsics * Σ_ext * transpose(jacobians.J_extrinsics)
    @test covariance ≈ transpose(covariance) atol=1e-12
    @test minimum(eigvals(Symmetric(covariance))) >= -1e-10
end

@testset "Implicit undistortion derivative" begin
    distortion = Float64[0.02, -0.003, 0.001, -0.002, 0.0004]
    xd, yd = 0.2, -0.14
    implicit_jacobian =
        UQ.undistortion_implicit_parameter_jacobian(xd, yd, distortion)
    finite_difference = Matrix{Float64}(undef, 2, 5)
    step = 1e-6
    for i in eachindex(distortion)
        plus = copy(distortion)
        minus = copy(distortion)
        plus[i] += step
        minus[i] -= step
        x_plus, y_plus = CAM.undistortion([xd], [yd], plus; max_iter=40, tol=1e-13)
        x_minus, y_minus = CAM.undistortion([xd], [yd], minus; max_iter=40, tol=1e-13)
        finite_difference[:, i] .=
            [x_plus[1] - x_minus[1], y_plus[1] - y_minus[1]] ./ (2 * step)
    end
    @test implicit_jacobian ≈ finite_difference atol=1e-7 rtol=2e-5
end

@testset "Relative-position assessment" begin
    J_offset = reshape(Float64[0.01, 0.0, 0.02, 0.0], 4, 1)
    normals = Float64[1 0; 1 0]
    joint_camera_normal_sigma = Float64[0.03, 0.06]
    t_mean = Float64[2.0, -1.0, 50.0]

    assessment = UQ.relative_position_assessment(
        J_offset, normals, joint_camera_normal_sigma, t_mean,
    )

    @test assessment.camera_marker_distance_mm == 50.0
    @test assessment.median_normal_sigma_per_mm_offset ≈ 0.015
    @test assessment.median_joint_camera_normal_sigma_mm ≈ 0.045
    @test assessment.break_even_offset_mm ≈ 3.0
    @test assessment.relative_scale_error_per_mm ≈ 1 / 50
    @test assessment.relative_surface_tension_error_per_mm ≈ 2 / 50
    @test_throws ArgumentError UQ.relative_position_assessment(
        J_offset, normals, joint_camera_normal_sigma, Float64[0.0, 0.0, 0.0],
    )
end

@testset "Shared-corner systematic sensitivity" begin
    intrinsics = Float64[800, 800, 320, 240, 0, 0, 0, 0, 0]
    K = Float64[800 0 320; 0 800 240; 0 0 1]
    corners = Float64[
        312 248
        328 248
        328 232
        312 232
    ]
    r = collect(range(-1.5, 1.5; length=21))
    z = 2 .- 0.6 .* r.^2
    virtual_u = 320 .+ 800 .* r ./ 500
    virtual_v = 240 .+ 800 .* z ./ 500
    assessment = UQ.shared_corner_shift_sensitivity(
        corners, K, zeros(5), 10.0, virtual_u, virtual_v, intrinsics,
        Matrix{Float64}(I, 3, 3), zeros(3), fill(0.006, length(r)),
    )

    @test length(assessment.cases) == 3
    @test assessment.median_joint_camera_sigma_um ≈ 6.0
    @test assessment.cases[1].translation_change_um[1] ≈ 62.5 atol=1e-2
    @test assessment.cases[1].translation_change_um[2] ≈ 0.0 atol=1e-2
    @test assessment.cases[2].translation_change_um[2] ≈ 62.5 atol=1e-2
    @test assessment.cases[2].translation_change_um[1] ≈ 0.0 atol=1e-2
    @test assessment.cases[3].translation_change_um[3] < 0
    @test all(
        isfinite,
        vcat(
            assessment.cases[1].rotation_change_mrad,
            assessment.cases[1].translation_change_um,
            assessment.cases[1].median_abs_normal_displacement_um,
            assessment.cases[1].apex_radius_change_um,
        ),
    )
end

const EDGE_TEST_ROOT = normpath(joinpath(@__DIR__, ".."))
const EDGE_TEST_IMAGE = joinpath(EDGE_TEST_ROOT, "data", "test_droplets", "zoom5.JPG")

function _with_edge_test_file(f::Function, data::Dict{String,Any})
    path, io = mktemp()
    close(io)
    try
        JLD2.jldopen(path, "w") do file
            for (key, value) in data
                file[key] = value
            end
        end
        return f(path)
    finally
        rm(path; force=true)
    end
end

function _valid_edge_test_data()
    return Dict{String,Any}(
        "u" => Float64[1000, 1010, 1020, 1030],
        "v" => Float64[1500, 1510, 1520, 1530],
        "sigma_n_px" => fill(0.3, 4),
        "sigma_description" => "1σ edge-fit uncertainty; estimated from fit residuals.",
        "image_file" => "data/test_droplets/zoom5.JPG",
        "student_note" => "extra keys are ignored",
    )
end

function _load_test_edge_data(path)
    return UQ.load_edge_uncertainty_data(
        path;
        image_root=EDGE_TEST_ROOT,
        expected_image_file=EDGE_TEST_IMAGE,
    )
end

@testset "Observed edge-data contract" begin
    valid = _valid_edge_test_data()
    _with_edge_test_file(valid) do path
        loaded = _load_test_edge_data(path)
        @test loaded.u == valid["u"]
        @test loaded.v == valid["v"]
        @test loaded.sigma_description == valid["sigma_description"]
        @test loaded.correlation_source == "uncorrelated"
        @test loaded.correlation == Matrix{Float64}(I, 4, 4)
    end

    mktempdir() do directory
        same_size_image = joinpath(directory, "drop_any_name.jpg")
        cp(EDGE_TEST_IMAGE, same_size_image)
        other_image = copy(valid)
        other_image["image_file"] = same_size_image
        _with_edge_test_file(other_image) do path
            @test _load_test_edge_data(path).image_path == realpath(same_size_image)
        end

        smaller_image = joinpath(directory, "small.png")
        FileIO.save(smaller_image, FileIO.load(EDGE_TEST_IMAGE)[1:10, 1:10])
        mismatched_dimensions = copy(valid)
        mismatched_dimensions["u"] = Float64[1, 2, 3, 4]
        mismatched_dimensions["v"] = Float64[1, 2, 3, 4]
        mismatched_dimensions["image_file"] = smaller_image
        _with_edge_test_file(mismatched_dimensions) do path
            @test_throws ArgumentError _load_test_edge_data(path)
        end
    end

    with_corr_length = copy(valid)
    with_corr_length["corr_length_px"] = 8.0
    _with_edge_test_file(with_corr_length) do path
        loaded = _load_test_edge_data(path)
        @test loaded.correlation_source == "corr_length_px"
        @test loaded.correlation[1, 2] ≈ exp(-hypot(10, 10) / 8)
    end

    with_corr_matrix = copy(valid)
    with_corr_matrix["corr_matrix"] = Matrix{Float64}(I, 4, 4)
    _with_edge_test_file(with_corr_matrix) do path
        @test _load_test_edge_data(path).correlation_source == "corr_matrix"
    end

    function rejects(data)
        _with_edge_test_file(data) do path
            @test_throws ArgumentError _load_test_edge_data(path)
        end
    end

    missing_key = copy(valid)
    delete!(missing_key, "sigma_n_px")
    _with_edge_test_file(missing_key) do path
        @test_throws ErrorException _load_test_edge_data(path)
    end
    missing_description = copy(valid)
    delete!(missing_description, "sigma_description")
    _with_edge_test_file(missing_description) do path
        @test_throws ErrorException _load_test_edge_data(path)
    end

    invalid = copy(valid)
    invalid["u"] = Float32[1000, 1010, 1020, 1030]
    rejects(invalid)
    invalid = copy(valid)
    invalid["v"] = Float64[1500, 1510, 1520]
    rejects(invalid)
    invalid = copy(valid)
    invalid["u"] = Float64[1000, NaN, 1020, 1030]
    rejects(invalid)
    invalid = copy(valid)
    invalid["u"] = Float64[0, 1010, 1020, 1030]
    rejects(invalid)
    invalid = copy(valid)
    invalid["sigma_n_px"] = Float64[0.3, 0.3, 0.0, 0.3]
    rejects(invalid)
    invalid = copy(valid)
    invalid["sigma_description"] = "  "
    rejects(invalid)
    invalid = copy(valid)
    invalid["sigma_description"] = 0.3
    rejects(invalid)
    invalid = copy(valid)
    invalid["corr_length_px"] = 0.0
    rejects(invalid)
    invalid = copy(valid)
    invalid["corr_length_px"] = 8.0
    invalid["corr_matrix"] = Matrix{Float64}(I, 4, 4)
    rejects(invalid)
    invalid = copy(valid)
    invalid["corr_matrix"] = Matrix{Float64}(I, 3, 3)
    rejects(invalid)
    invalid = copy(valid)
    invalid["corr_matrix"] = [1.0 0.2 0 0; 0.1 1 0 0; 0 0 1 0; 0 0 0 1]
    rejects(invalid)
    invalid = copy(valid)
    invalid["corr_matrix"] = [1.0 0 0 0; 0 0.9 0 0; 0 0 1 0; 0 0 0 1]
    rejects(invalid)
    invalid = copy(valid)
    invalid["corr_matrix"] = [1.0 2 0 0; 2 1 0 0; 0 0 1 0; 0 0 0 1]
    rejects(invalid)

    accepted_r = Float64[0.01, 0.25, 0.7, 1.4, 0.915]
    accepted_z = Float64[-4.0, -3.0, -2.0, -1.0, 0.0]
    @test UQ.validate_observed_contour_order(accepted_r, accepted_z, 0.915) === nothing
    @test_throws ArgumentError UQ.validate_observed_contour_order(
        reverse(accepted_r), reverse(accepted_z), 0.915,
    )
    @test_throws ArgumentError UQ.validate_observed_contour_order(
        Float64[0.01, -0.4, -0.8, -0.3, 0.915], accepted_z, 0.915,
    )
    @test UQ.validate_observed_contour_order(
        -accepted_r, accepted_z, 0.915,
    ) === nothing
end

@testset "Observed edge covariance plumbing" begin
    intrinsics = Float64[800, 800, 2304, 1728, 0, 0, 0, 0, 0]
    R = Matrix{Float64}(I, 3, 3)
    t = Float64[0, 0, 500]
    t_rel = zeros(3)
    r = collect(range(0, 0.915; length=5))
    z = collect(range(-4, 0; length=5))
    parameters = vcat(intrinsics, zeros(6), 0.0)
    virtual_u, virtual_v, _ = UQ.forward_project_camtodrop(
        r, z, parameters, R, t, R, t_rel;
        plane_coordinate=0.0,
    )
    edge_data = _valid_edge_test_data()
    edge_data["u"] = virtual_u .+ 1
    edge_data["v"] = virtual_v .+ 1
    edge_data["sigma_n_px"] = fill(0.3, length(r))

    _with_edge_test_file(edge_data) do path
        loaded = _load_test_edge_data(path)
        @test loaded.u .- 1 ≈ virtual_u
        @test loaded.v .- 1 ≈ virtual_v

        jacobian_virtual = UQ.camtodrop_jacobians(
            virtual_u, virtual_v, intrinsics, R, t, R, t_rel;
            plane_coordinate=0.0,
        )
        jacobian_observed = UQ.camtodrop_jacobians(
            loaded.u .- 1, loaded.v .- 1, intrinsics, R, t, R, t_rel;
            plane_coordinate=0.0,
        )
        Σ_intrinsics = Matrix(Diagonal(fill(1e-8, 9)))
        Σ_ext = Matrix(Diagonal(fill(1e-6, 6)))
        terms_virtual = UQ.covariance_terms(
            jacobian_virtual, Σ_intrinsics, Σ_ext, Σ_ext, 0.3, 0.0,
        )
        terms_observed = UQ.covariance_terms(
            jacobian_observed, Σ_intrinsics, Σ_ext, Σ_ext, 0.3, 0.0,
        )
        @test terms_observed.camera_full ≈ terms_virtual.camera_full

        normals = UQ.contour_normals(
            jacobian_observed.coordinates[1:2:end],
            jacobian_observed.coordinates[2:2:end],
        )
        edge_pixel_covariance = UQ.edge_normal_pixel_covariance(
            loaded.u, loaded.v, loaded.sigma_n_px, loaded.correlation,
        )
        edge_model_covariance =
            jacobian_observed.J_pixels * edge_pixel_covariance *
            transpose(jacobian_observed.J_pixels)
        Σ_edge_normal =
            UQ.normal_covariance(edge_model_covariance, normals)
        Σ_isotropic_normal = UQ.normal_covariance(terms_observed.pixel, normals)
        @test Σ_edge_normal ≈ Σ_isotropic_normal atol=1e-10 rtol=1e-8
    end

    R_mirror = Matrix(Diagonal(Float64[-1, 1, -1]))
    virtual_left_u, _, _ = UQ.forward_project_camtodrop(
        r, z, parameters, R_mirror, t, R, t_rel;
        plane_coordinate=0.0,
    )
    mirrored_u, mirrored_v, _ = UQ.forward_project_camtodrop(
        -r, z, parameters, R_mirror, t, R, t_rel;
        plane_coordinate=0.0,
    )
    @test all(virtual_left_u[2:end] .< intrinsics[3])
    @test all(mirrored_u[2:end] .> intrinsics[3])
    mirrored_edge_data = _valid_edge_test_data()
    mirrored_edge_data["u"] = mirrored_u .+ 1
    mirrored_edge_data["v"] = mirrored_v .+ 1
    mirrored_edge_data["sigma_n_px"] = fill(0.3, length(r))
    _with_edge_test_file(mirrored_edge_data) do path
        loaded = _load_test_edge_data(path)
        recovered_mirrored_r, recovered_mirrored_z = UQ.inverse_project_camtodrop(
            loaded.u .- 1, loaded.v .- 1, parameters, R_mirror, t, R, t_rel;
            plane_coordinate=0.0,
        )
        @test recovered_mirrored_r ≈ -r atol=1e-8
        @test recovered_mirrored_z ≈ z atol=1e-8
        @test UQ.validate_observed_contour_order(
            recovered_mirrored_r, recovered_mirrored_z, 0.915,
        ) === nothing
        mirrored_jacobians = UQ.camtodrop_jacobians(
            loaded.u .- 1, loaded.v .- 1, intrinsics, R_mirror, t, R, t_rel;
            plane_coordinate=0.0,
        )
        mirrored_normals = UQ.contour_normals(
            mirrored_jacobians.coordinates[1:2:end],
            mirrored_jacobians.coordinates[2:2:end],
        )
        @test all(isfinite, mirrored_jacobians.J_pixels)
        @test all(isfinite, UQ.normal_covariance(
            mirrored_jacobians.J_pixels *
            Matrix{Float64}(I, size(mirrored_jacobians.J_pixels, 2),
                            size(mirrored_jacobians.J_pixels, 2)) *
            transpose(mirrored_jacobians.J_pixels),
            mirrored_normals,
        ))
    end
end
