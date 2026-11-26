using ForwardDiff
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
