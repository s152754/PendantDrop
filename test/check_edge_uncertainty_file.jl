using JLD2
using LinearAlgebra

include("../src/forwardmodel.jl")

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const CAMERA = Model.Camera
const UQ = CAMERA.IntrinsicUncertainty
const EDGE_DATA_FILE = isempty(ARGS) ?
    joinpath(PROJECT_ROOT, "docs", "example_edge_data.jld2") :
    abspath(ARGS[1])
const INTRINSIC_FILE = joinpath(
    PROJECT_ROOT, "data", "calibration_camera", "zoom5", "detected",
    "camint_all.jld2",
)
const POSE_RESULTS_FILE = joinpath(
    PROJECT_ROOT, "data", "test_droplets", "camtodrop_uncertainty",
    "uncertainty_camtodrop.jld2",
)
const IMAGE_FILE = joinpath(PROJECT_ROOT, "data", "test_droplets", "zoom5.JPG")

edge = UQ.load_edge_uncertainty_data(
    EDGE_DATA_FILE;
    image_root=PROJECT_ROOT,
    expected_image_file=IMAGE_FILE,
)
intrinsics_data = JLD2.load(INTRINSIC_FILE)
pose_data = JLD2.load(POSE_RESULTS_FILE)
UQ.require_keys(intrinsics_data, ("camera_matrix", "dist_coeffs"), INTRINSIC_FILE)
UQ.require_keys(pose_data, ("rvec_mean", "t_mean"), POSE_RESULTS_FILE)

K = Matrix{Float64}(intrinsics_data["camera_matrix"])
dist = vec(Float64.(intrinsics_data["dist_coeffs"]))[1:5]
intrinsics = Float64[
    K[1, 1], K[2, 2], K[1, 3], K[2, 3], dist[1], dist[2],
    dist[3], dist[4], dist[5],
]
R_mean = UQ._right_perturbed_rotation(
    Matrix{Float64}(I, 3, 3), Float64.(pose_data["rvec_mean"]),
)
r, z = UQ.inverse_project_camtodrop(
    edge.u .- 1, edge.v .- 1, vcat(intrinsics, zeros(6), 0.0),
    R_mean, Float64.(pose_data["t_mean"]),
    CAMERA.getRrel_fn(0.0, 0.0, 0.0), Float64[-8.3, 0.6, 0.0];
    plane_coordinate=0.0,
)
UQ.validate_observed_contour_order(r, z, 1.83 / 2)

println("Sigma description (1σ): ", edge.sigma_description)
println(
    "Edge-data contract, image dimensions, and apex-to-needle ordering passed: $EDGE_DATA_FILE",
)
