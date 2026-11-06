"""
Camera calibration utilities

This include-only file provides checkerboard image selection and calibration,
intrinsic-parameter uncertainty estimation, and plotting helpers. It relies on
packages and helpers loaded by the parent camera module; do not include it alone.
"""


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
function detect_aruco_corners(img_file::String, K, dist;
                              marker_dict=nothing, aruco_params=nothing)
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