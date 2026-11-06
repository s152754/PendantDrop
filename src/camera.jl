module Camera

using JLD2
using PyCall
using FileIO
using Images
using Plots
using LinearAlgebra
using Statistics
using Random

export ArUcopose, calibrate_camera, summarize_camera_calibration_results, uncertainty_camera_intrinsics, uncertainty_camera_intrinsics_by_subset, calibrate_camera_from_images, random_image_subset, determine_random_image_subsets, list_camera_images, project_model_to_image, pixel_to_model, getRrel_fn, R_to_euler_zyx
export IntrinsicUncertainty

function Forward_projection(rm, zm, R, t, K, dist_coeffs, R_rel, t_rel)

    u, v, - = Camera.project_model_to_image(rm, zm, sm, R, t, K, dist_coeffs, R_rel, t_rel)

    return u, v

end

function Inverse_projection(u, v, R, t, K, dist_coeffs, R_rel, t_rel)

    return rm, zm

end

function calibrate_camera(img_dir::String, out_dir::String; checkerboard_size::Tuple{Int,Int}=(13, 9), square_size::Float64=0.021, img_files::Union{Nothing, Vector{String}} = nothing)
    """
    Calibrate camera from checkerboard images.
    
    Args:
        img_dir:            Directory containing checkerboard calibration images
        out_dir:            Directory to save detected images and calibration files
        checkerboard_size:  Tuple of (horizontal_corners, vertical_corners)             [-]
        square_size:        Physical size of checkerboard squares                       [m]
    
    Returns:

        Dict with calibration results (K, distortion coefficients and rms) or nothing if calibration fails
    """
    ensure_directory(out_dir)

    result = py"""
    import cv2, numpy as np, glob, os, sys

    # --- Configuration ---
    checker = tuple($checkerboard_size)
    square  = $square_size
    os.makedirs($out_dir, exist_ok=True)

    # --- Load all images ---
    # img_files = glob.glob(os.path.join($img_dir, "*.JPG")) + glob.glob(os.path.join($img_dir, "*.jpg"))

    # --- Load images ---

    if $img_files is None:
        img_files = (glob.glob(os.path.join($img_dir, "*.JPG")) + glob.glob(os.path.join($img_dir, "*.jpg")))
        print(f"Using all images in directory: {len(img_files)} found")

    else:
        img_files = list($img_files)
        print(f"Using supplied image list: {len(img_files)} files")

    # --- Initialize storage for corner points ---
    objpoints   = []  # 3D points in real world space
    imgpoints   = []  # 2D points in image plane
    saved_paths = []

    # --- Create template object points for checkerboard ---
    objp        = np.zeros((checker[0] * checker[1], 3), np.float32)
    objp[:, :2] = np.mgrid[0:checker[0], 0:checker[1]].T.reshape(-1, 2)
    objp        *= square

    # --- Detect corners in each image ---
    for img_path in img_files:
        img = cv2.imread(img_path)
        if img is None:
            print(f"Could not read {img_path}")
            sys.stdout.flush()
            continue

        # --- Normalize to 3-channel BGR ---
        if img.ndim == 2:
            img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
        elif img.shape[2] == 4:
            img = cv2.cvtColor(img, cv2.COLOR_RGBA2BGR)

        # --- Convert to grayscale ---
        gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
        
        # --- Find checkerboard corners ---
        flags           = cv2.CALIB_CB_ADAPTIVE_THRESH + cv2.CALIB_CB_NORMALIZE_IMAGE + cv2.CALIB_CB_FAST_CHECK
        # ret, corners    = cv2.findChessboardCorners(gray, checker, flags=flags)
        ret, corners = cv2.findChessboardCornersSB(gray, checker)
        print(f"  corners detected: {ret}")

        # --- Refine corners if detected ---
        if ret:
            corners         = corners.astype(np.float32)
            criteria        = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 30, 0.001)
            corners_refined = cv2.cornerSubPix(gray, corners, (11, 11), (-1, -1), criteria)

            objpoints.append(objp)
            imgpoints.append(corners_refined)

            # --- Draw and save image with detected corners ---
            cv2.drawChessboardCorners(img, checker, corners_refined, ret)
            out_img_path = os.path.join($out_dir, os.path.basename(img_path))
            cv2.imwrite(out_img_path, img)
            saved_paths.append(out_img_path)
        else:
            print(f"No corners in {os.path.basename(img_path)}")
            sys.stdout.flush()

        print(f"Detected corners in {len(saved_paths)} images")
        print(f"objpoints count: {len(objpoints)}, imgpoints count: {len(imgpoints)}")

    # --- Calibrate camera if corners were found ---
    if len(objpoints) > 0:
        h, w = gray.shape[:2]
        print(f"Calibrating with {len(objpoints)} valid views")
        ret, mtx, dist, rvecs, tvecs = cv2.calibrateCamera(objpoints, imgpoints, (w, h), None, None)

        # --- RMS reprojection error per image ---
        image_rms = []

        for i in range(len(objpoints)):

            projected, _ = cv2.projectPoints(
                objpoints[i],
                rvecs[i],
                tvecs[i],
                mtx,
                dist
            )

            residuals = (
                imgpoints[i].reshape(-1, 2)
                -
                projected.reshape(-1, 2)
            )

            rms_i = np.sqrt(
                np.mean(
                    residuals[:, 0]**2 +
                    residuals[:, 1]**2
                )
            )

            image_rms.append(float(rms_i))

            print(
                f"Image {i+1:02d}: RMS reprojection error = "
                f"{rms_i:.3f} px"
            )

        # # --- Save calibration matrices ---
        # np.save(os.path.join($out_dir, "camera_matrix.npy"), mtx)
        # np.save(os.path.join($out_dir, "dist_coeffs.npy"), dist)

        print("Calibration RMS:", ret)
        print("Camera matrix:\n", mtx)
        print("Distortion coeffs:\n", dist)

        rms = float(ret)
        cammat = mtx.tolist()
        distcoef = np.array(dist).reshape(-1).tolist()
        detpath = list(saved_paths)

        result = (rms, cammat, distcoef, detpath, 
        [o.tolist() for o in objpoints],
        [i.reshape(-1,2).tolist() for i in imgpoints],
        [r.reshape(-1).tolist() for r in rvecs],
        [t.reshape(-1).tolist() for t in tvecs],
        (w, h),
        image_rms)
    else:
        print("No corners detected; calibration skipped")
        result = None
        print(f"werkt result = {result}")

    result

    """

    return py"result"
end

# inverse projection

function undistortion(xd, yd, dist_coeffs;
                           max_iter=10,
                           tol=1e-8())

    k1, k2, p1, p2, k3 = dist_coeffs

    # Initial guess
    x = copy(xd)
    y = copy(yd)

    for _ in 1:max_iter

        r2 = x.^2 .+ y.^2
        r4 = r2.^2
        r6 = r2.^3

        radial = 1 .+ k1 .* r2 .+ k2 .* r4 .+ k3 .* r6

        Δx = 2 .* p1 .* x .* y .+
             p2 .* (r2 .+ 2 .* x.^2)

        Δy = p1 .* (r2 .+ 2 .* y.^2) .+
             2 .* p2 .* x .* y

        x_new = (xd .- Δx) ./ radial
        y_new = (yd .- Δy) ./ radial

        err = maximum(sqrt.((x_new .- x).^2 .+ (y_new .- y).^2))

        x .= x_new
        y .= y_new

        if err < tol
            break
        end
    end

    return x, y

end

# forward projection
 function distortion(x, y, dist_coeffs::AbstractVector)
    """
    Apply distortion to normalized image coordinates (x, y) using distortion coefficients.

    Arg:
    - x, y         : normalized image coordinates
    - dist_coeffs  : distortion coefficients [k1, k2, p1, p2, k3]

    return:
    - x_dist, y_dist : distorted normalized image coordinates
    """
    k1, k2, p1, p2, k3 = dist_coeffs
    r2 = x.^2 .+ y.^2
    r4 = r2.^2
    r6 = r2.^3

    radial = 1 .+ k1*r2 .+ k2*r4 .+ k3*r6

    x_radial = x .* radial
    y_radial = y .* radial

    x_tangential = 2p1 .* x .* y .+ p2 .* (r2 .+ 2 .* x.^2)
    y_tangential = p1 .* (r2 .+ 2 .* y.^2) .+ 2p2 .* x .* y

    x_dist = x_radial .+ x_tangential
    y_dist = y_radial .+ y_tangential

    
    return x_dist, y_dist
end


function ArUcopose(K, dist_coeffs, img_file, marker_length_mm; return_plot::Bool=false, img_dir=nothing)
    """
    Detect ArUco markers in a single image and return pose information.

    Arguments:
    - K         : camera matrix
    - dist_coeffs : distortion coefficients
    - img_file  : path to a single image file

    Keyword arguments:
    - marker_length_mm: physical marker size in millimeters
    - img_dir: optional directory to prepend if img_file is a relative filename

    Returns: (R, tvec, detection_result)
    """

    img_path = string(img_file)
    if !isabspath(img_path)
        candidates = [joinpath(pwd(), img_path)]
        if img_dir !== nothing && !isempty(img_dir)
            push!(candidates, joinpath(img_dir, img_path))
        end
        push!(candidates, joinpath(pwd(), "test_images", "detected", img_path))
        for candidate in candidates
            if isfile(candidate)
                img_path = candidate
                break
            end
        end
    end

    py"""
    import cv2, numpy as np, os, sys

    K                   = np.array($K)
    dist                = np.array($dist_coeffs).flatten()
    marker_length_mm    = $marker_length_mm
    img_path            = $img_path

    if not os.path.exists(img_path):
        raise FileNotFoundError(f"Image file not found: {img_path}")

    img = cv2.imread(img_path)
    if img is None:
        raise ValueError(f"Could not read {img_path}")

    if img.ndim == 2:
        img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
    elif img.shape[2] == 4:
        img = cv2.cvtColor(img, cv2.COLOR_RGBA2BGR)

    gray        = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    aruco_dict  = cv2.aruco.getPredefinedDictionary(cv2.aruco.DICT_4X4_50)
    params      = cv2.aruco.DetectorParameters()

    ### added by Ariana
    params.adaptiveThreshWinSizeMin  = 15
    params.adaptiveThreshWinSizeMax  = 99     # push this well above one module's pixel size
    params.adaptiveThreshWinSizeStep = 10     # coarser steps, fewer windows
    params.minMarkerPerimeterRate = 0.05      # less noise
    ###

    detector    = cv2.aruco.ArucoDetector(aruco_dict, params)

    corners, ids, rejected = detector.detectMarkers(gray)
    if ids is None:
        raise ValueError(f"No markers detected in {img_path}")

    obj_points = np.array([
        [-marker_length_mm/2,  marker_length_mm/2, 0],
        [ marker_length_mm/2,  marker_length_mm/2, 0],
        [ marker_length_mm/2, -marker_length_mm/2, 0],
        [-marker_length_mm/2, -marker_length_mm/2, 0]
    ], dtype=np.float32)

    img_points = corners[0].reshape(4, 2).astype(np.float32)
    success, rvec, tvec = cv2.solvePnP(obj_points, img_points, K, dist)
    if not success:
        raise RuntimeError("solvePnP failed")

    R, _ = cv2.Rodrigues(rvec)
    marker_dist = float(tvec[2][0])

    proj_points, _ = cv2.projectPoints(obj_points, rvec, tvec, K, dist)
    proj_points = proj_points.reshape(-1, 2)
    px_dist = np.linalg.norm(proj_points[0] - proj_points[1])
    px_per_mm = px_dist / marker_length_mm
    mm_per_px = marker_length_mm / px_dist

    detection_result = {
        "filename": os.path.basename(img_path),
        "marker_ids": ids.flatten().tolist(),
        "marker_dist_mm": marker_dist,
        "px_per_mm": float(px_per_mm),
        "mm_per_px": float(mm_per_px),
    }

    # (R, tvec, detection_result)
    """

    # Try to create the grayscale plot (useful for quick visual check)
    p = nothing
    try
        # load via Julia and convert to Float grayscale
        jl_img = load(img_path)
        jl_gray = Float32.(Gray.(jl_img))
        p = heatmap(reverse(eachcol(jl_gray)); color=:gray, axis=false, legend=false, title="Gray: $(basename(img_path))")
    catch e
        @warn "Could not create grayscale plot: $e"
    end

    if return_plot
        return Array(py"R"), vec(Array(py"tvec")), py"detection_result", p
    end

    return Array(py"R"), vec(Array(py"tvec")), py"detection_result"
end


function ensure_directory(path::String)
    if !isdir(path)
        mkpath(path)
    end
end

function _coerce_camera_matrix(Kmat)
    if Kmat === nothing
        return zeros(3, 3)
    elseif Kmat isa AbstractMatrix
        return Float64.(Matrix(Kmat))
    elseif Kmat isa AbstractVector
        if length(Kmat) == 9
            return reshape(Float64.(Kmat), 3, 3)
        elseif !isempty(Kmat) && all(x -> x isa AbstractVector || x isa AbstractMatrix, Kmat)
            rows = [reshape(Float64.(collect(row)), 1, length(row)) for row in Kmat]
            return reduce(vcat, rows)
        end
    end

    return Float64.(Matrix(Kmat))
end

function _coerce_distortion_coeffs(distcoef)
    if distcoef === nothing
        return Float64[]
    elseif distcoef isa AbstractVector
        if !isempty(distcoef) && distcoef[1] isa AbstractVector
            return Float64.(collect(Iterators.flatten(distcoef)))
        end
        return Float64.(distcoef)
    elseif distcoef isa AbstractMatrix
        return vec(Float64.(distcoef))
    end

    return Float64[Float64(distcoef)]
end

function summarize_camera_calibration_results(calibration_results)
    """
    Summarize intrinsic calibration parameters from a collection of calibration runs.

    The input can contain any mix of:
    - tuples returned by `calibrate_camera`: `(rms, Kmat, distcoef, images)`
    - Dicts with `"camera_matrix"`, `"dist_coeffs"`, and `"rms"`
    - `nothing` or other invalid entries, which are ignored
    """

    if calibration_results isa Dict || calibration_results isa Tuple
        calibration_results = [calibration_results]
    end

    if isempty(calibration_results)
        return Dict(
            "K_mean" => zeros(3, 3),
            "K_std" => zeros(3, 3),
            "dist_mean" => Float64[],
            "dist_std" => Float64[],
            "rms_mean" => NaN,
            "rms_std" => NaN,
            "n_results" => 0,
        )
    end

    K_samples = Matrix{Float64}[]
    dist_samples = Matrix{Float64}[]
    rms_values = Float64[]

    for result in calibration_results
        if result === nothing
            continue
        end

        try
            if result isa Dict
                Kmat = get(result, "camera_matrix", get(result, "K", nothing))
                distcoef = get(result, "dist_coeffs", get(result, "dist", nothing))
                rms = get(result, "rms", nothing)
            elseif result isa Tuple && length(result) >= 3
                rms = result[1]
                Kmat = result[2]
                distcoef = result[3]
            else
                continue
            end

            if Kmat === nothing
                continue
            end

            push!(K_samples, _coerce_camera_matrix(Kmat))
            push!(dist_samples, _coerce_distortion_coeffs(distcoef))

            if rms !== nothing
                push!(rms_values, Float64(rms))
            end
        catch
            continue
        end
    end

    if isempty(K_samples)
        K_mean = zeros(3, 3)
        K_std = zeros(3, 3)
    else
        K_mean = zeros(3, 3)
        K_std = zeros(3, 3)
        for i in 1:3, j in 1:3
            vals = [K[i, j] for K in K_samples]
            K_mean[i, j] = mean(vals)
            K_std[i, j] = std(vals)
        end
    end

    if isempty(dist_samples)
        dist_mean = Float64[]
        dist_std = Float64[]
    else
        n_coeffs = length(dist_samples[1])
        dist_mean = zeros(n_coeffs)
        dist_std = zeros(n_coeffs)
        for k in 1:n_coeffs
            vals = [d[k] for d in dist_samples]
            dist_mean[k] = mean(vals)
            dist_std[k] = std(vals)
        end
    end

    return Dict(
        "K_mean" => K_mean,
        "K_std" => K_std,
        "dist_mean" => dist_mean,
        "dist_std" => dist_std,
        "rms_mean" => isempty(rms_values) ? NaN : mean(rms_values),
        "rms_std" => isempty(rms_values) ? NaN : std(rms_values),
        "n_results" => length(K_samples),
    )
end



function list_camera_images(img_dir::String; image_extensions::Vector{String}=String[".JPG"])::Vector{String}
    if !isdir(img_dir)
        return String[]
    end

    image_files = String[]
    for name in readdir(img_dir)
        path = joinpath(img_dir, name)
        if isfile(path)
            lower_name = lowercase(name)
            if any(endswith(lower_name, lowercase(ext)) for ext in image_extensions)
                push!(image_files, abspath(path))
            end
        end
    end

    return sort(image_files)
end

function random_image_subset(img_dir::String; subset_size::Int=5, seed::Union{Int,Nothing}=nothing, image_extensions::Vector{String}=String[".JPG"])::Vector{String}
    image_files = list_camera_images(img_dir; image_extensions=image_extensions)
    if isempty(image_files)
        return String[]
    end

    if subset_size <= 0
        return String[]
    end

    shuffled = copy(image_files)
    if seed !== nothing
        rng = MersenneTwister(seed)
        Random.shuffle!(rng, shuffled)
    else
        Random.shuffle!(shuffled)
    end

    return shuffled[1:min(subset_size, length(shuffled))]
end

function determine_random_image_subsets(img_dir::String; n_images_per_subset::Int=5, seed::Union{Int,Nothing}=nothing, image_extensions::Vector{String}=String[".JPG"])
    image_files = list_camera_images(img_dir; image_extensions=image_extensions)
    if isempty(image_files)
        return Vector{Vector{String}}()
    end

    shuffled = copy(image_files)
    if seed !== nothing
        rng = MersenneTwister(seed)
        Random.shuffle!(rng, shuffled)
    else
        Random.shuffle!(shuffled)
    end

    subsets = Vector{Vector{String}}()
    for start in 1:n_images_per_subset:length(shuffled)
        stop = min(start + n_images_per_subset - 1, length(shuffled))
        push!(subsets, shuffled[start:stop])
    end

    return subsets
end

function project_checkerboard_points(
    objpoints,
    rvecs,
    tvecs,
    K,
    dist
)

    py"""
import cv2
import numpy as np

all_projected = []

for objp, rvec, tvec in zip(
        $objpoints,
        $rvecs,
        $tvecs):

    objp = np.array(objp, dtype=np.float32)

    rvec = np.array(rvec, dtype=np.float64).reshape(3,1)
    tvec = np.array(tvec, dtype=np.float64).reshape(3,1)

    projected, _ = cv2.projectPoints(
        objp,
        rvec,
        tvec,
        np.array($K, dtype=np.float64),
        np.array($dist, dtype=np.float64)
    )

    all_projected.append(
        projected.reshape(-1,2).tolist()
    )

result = all_projected
"""

    return py"result"

end


function compute_reprojection_residuals(
    objpoints,
    imgpoints,
    rvecs,
    tvecs,
    K,
    dist
)

    py"""
import cv2
import numpy as np

all_residuals = []

for objp, imgp, rvec, tvec in zip(
        $objpoints,
        $imgpoints,
        $rvecs,
        $tvecs):

    objp = np.array(objp, dtype=np.float32)
    imgp = np.array(imgp, dtype=np.float64)

    rvec = np.array(rvec, dtype=np.float64).reshape(3,1)
    tvec = np.array(tvec, dtype=np.float64).reshape(3,1)

    projected, _ = cv2.projectPoints(
        objp,
        rvec,
        tvec,
        np.array($K, dtype=np.float64),
        np.array($dist, dtype=np.float64)
    )

    projected = projected.reshape(-1,2)

    residuals = imgp - projected

    all_residuals.append(residuals)

all_residuals = np.vstack(all_residuals)

sigma_x = float(np.std(all_residuals[:,0]))
sigma_y = float(np.std(all_residuals[:,1]))

result = (
    sigma_x,
    sigma_y,
    all_residuals.tolist()
)
"""
    
    return py"result"

end

function bootstrap_calibration(
    objpoints,
    imgpoints,
    image_size,
    sigma_x,
    sigma_y;
    N_boot = 500
)

    py"""
import cv2
import numpy as np

objpoints  = $objpoints
imgpoints  = $imgpoints
image_size = tuple($image_size)

sigma_x = float($sigma_x)
sigma_y = float($sigma_y)

N_boot = int($N_boot)

samples = []

for boot in range(N_boot):

    perturbed_imgpoints = []

    # Perturb detected corners
    for corners in imgpoints:

        corners = np.array(
            corners,
            dtype=np.float32
        ).reshape(-1,2)

        noise = np.zeros_like(corners)

        noise[:,0] = np.random.normal(
            0,
            sigma_x,
            size=len(corners)
        )

        noise[:,1] = np.random.normal(
            0,
            sigma_y,
            size=len(corners)
        )

        corners_perturbed = (
            corners + noise
        ).astype(np.float32)

        perturbed_imgpoints.append(
            corners_perturbed.reshape(-1,1,2)
        )

    objpoints_cv = [
        np.array(o, dtype=np.float32)
        for o in objpoints
    ]

    ret, K, dist, rvecs, tvecs = cv2.calibrateCamera(
        objpoints_cv,
        perturbed_imgpoints,
        image_size,
        None,
        None
    )

    dist = np.array(dist).reshape(-1)

    fx = K[0,0]
    fy = K[1,1]
    cx = K[0,2]
    cy = K[1,2]

    phi = [fx, fy, cx, cy]
    phi.extend(dist.tolist())

    samples.append(phi)

samples = np.array(samples)

mean_phi = np.mean(
    samples,
    axis=0
)

cov_phi = np.cov(
    samples,
    rowvar=False
)

result = (
    samples.tolist(),
    mean_phi.tolist(),
    cov_phi.tolist()
)
"""
    
    return py"result"

end


function calibrate_camera_from_images(image_files::Vector{String}; checkerboard_size::Tuple{Int,Int}=(13, 9), square_size::Float64=0.021, out_dir::Union{String,Nothing}=nothing, reference_K=nothing, reference_dist=nothing)
    if isempty(image_files)
        return Dict(
            "camera_matrix" => nothing,
            "dist_coeffs" => nothing,
            "rms" => nothing,
            "images" => String[],
            "n_images" => 0,
            "success" => false,
        )
    end

    if out_dir !== nothing
        ensure_directory(out_dir)
    end

    result = nothing
    try
        result = py"""
        import cv2, numpy as np, os, sys

        checker = tuple($checkerboard_size)
        square = float($square_size)
        image_files = $image_files
        out_dir = None if $out_dir is None else $out_dir

        if out_dir is not None:
            os.makedirs(out_dir, exist_ok=True)

        objpoints = []
        imgpoints = []
        used_images = []
        saved_paths = []

        objp = np.zeros((checker[0] * checker[1], 3), np.float32)
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

            gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
            last_gray = gray

            ret, corners = cv2.findChessboardCornersSB(gray, checker)
            if not ret:
                continue

            corners = corners.astype(np.float32)
            criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 30, 0.001)
            corners_refined = cv2.cornerSubPix(gray, corners, (11, 11), (-1, -1), criteria)

            objpoints.append(objp)
            imgpoints.append(corners_refined)
            used_images.append(os.path.basename(img_path))

            if out_dir is not None:
                out_img_path = os.path.join(out_dir, os.path.basename(img_path))
                cv2.drawChessboardCorners(img, checker, corners_refined, True)
                cv2.imwrite(out_img_path, img)
                saved_paths.append(out_img_path)

        if len(objpoints) > 0 and last_gray is not None:
            h, w = last_gray.shape[:2]
            init_K = None if $reference_K is None else np.array($reference_K, dtype=np.float64)
            init_dist = None if $reference_dist is None else np.array($reference_dist, dtype=np.float64).reshape(1, -1)
            flags = cv2.CALIB_USE_INTRINSIC_GUESS if init_K is not None else 0
            ret, mtx, dist, _, _ = cv2.calibrateCamera(objpoints, imgpoints, (w, h), init_K, init_dist, flags=flags)
            result = (float(ret), mtx.tolist(), dist.reshape(-1).tolist(), used_images, len(used_images), True)
        else:
            result = (None, None, None, used_images, len(used_images), False)

        result
        """
        result = py"result"
    catch err
        @warn "Calibration failed for the selected image subset: $err"
        result = nothing
    end

    if result === nothing
        return Dict(
            "camera_matrix" => nothing,
            "dist_coeffs" => nothing,
            "rms" => nothing,
            "images" => String[],
            "n_images" => 0,
            "success" => false,
        )
    end

    if result isa Tuple && length(result) == 6
        rms, cammat, distcoef, used_images, n_images, success = result
    else
        return Dict(
            "camera_matrix" => nothing,
            "dist_coeffs" => nothing,
            "rms" => nothing,
            "images" => String[],
            "n_images" => 0,
            "success" => false,
        )
    end

    return Dict(
        "camera_matrix" => cammat,
        "dist_coeffs" => distcoef,
        "rms" => rms,
        "images" => used_images,
        "n_images" => n_images,
        "success" => success,
    )
end

function project_image_to_model()

    return (rm , zm , sm)
end 

function project_model_to_image(rm::Vector{Float64}, zm::Vector{Float64}, sm::Vector{Float64}, R::Matrix{Float64}, t::Vector{Float64}, K::Matrix{Float64}, dist_coeffs::AbstractVector, R_rel::Matrix{Float64}, t_rel::Vector{Float64})
    """
    Project 3D model points to image space.

    Arg:
    - rm, zm, sm   : model coordinates (radial, axial, arclength)
    - R, t         : Rotation matrix (R) and translation vector (t) for pose of ArUco marker in camera frame, 
    - R_rel, t_rel : Rotation matrix (R_rel) and translation vector (t_rel) for pose of object in marker frame
    - X            : 4×N homogeneous object points (r, z, s, 1) [mm?]
    - K            : camera intrinsic matrix
    - dist_coeffs  : distortion coefficients

    return:
    - u, v, w      : pixel coordinates (u, v) and depth coordinate (w) along the camera optical axis
    """

    @assert length(rm) == length(zm) == length(sm) "rm, zm, and sm must have same length"
    n = length(rm)
    X = vcat(reshape(rm, 1, n), reshape(zm, 1, n), reshape(sm, 1, n), ones(1, n))

    # Compose transform: object → camera
    R_total = R * R_rel
    t_total = R * t_rel + t

    # Transform points into camera frame
    X_cam = R_total * X[1:3, :] .+ t_total

    # Normalize
    x = X_cam[1, :] ./ X_cam[3, :]
    y = X_cam[2, :] ./ X_cam[3, :]

    # Distortion
    x_dist, y_dist = distortion(x, y, dist_coeffs)

    # Apply intrinsics
    fx = K[1,1]
    fy = K[2,2]
    cx = K[1,3]
    cy = K[2,3]

    u = fx .* x_dist .+ cx
    v = fy .* y_dist .+ cy
    w = X_cam[3, :]

    if length(u) == 1
        return (u[1], v[1], w[1])
    end

    return (u, v, w)
end

function pixel_to_model(
    u::AbstractVector,
    v::AbstractVector,
    R::AbstractMatrix,
    t::AbstractVector,
    K::AbstractMatrix,
    dist_coeffs::AbstractVector,
    R_rel::AbstractMatrix,
    t_rel::AbstractVector;
    plane_coordinate::Real = 0.0,
    iterations::Int = 8,
)
    length(u) == length(v) ||
        throw(ArgumentError("u and v must have the same length"))
    iterations >= 1 || throw(ArgumentError("iterations must be positive"))

    fx, fy = Float64(K[1, 1]), Float64(K[2, 2])
    cx, cy = Float64(K[1, 3]), Float64(K[2, 3])
    xd = (Float64.(u) .- cx) ./ fx
    yd = (Float64.(v) .- cy) ./ fy

    k1, k2, p1, p2, k3 = Float64.(dist_coeffs)
    x = copy(xd)
    y = copy(yd)
    for _ in 1:iterations
        r2 = x.^2 .+ y.^2
        radial = 1 .+ k1 .* r2 .+ k2 .* r2.^2 .+ k3 .* r2.^3
        x_tangential = 2 .* p1 .* x .* y .+ p2 .* (r2 .+ 2 .* x.^2)
        y_tangential = p1 .* (r2 .+ 2 .* y.^2) .+ 2 .* p2 .* x .* y
        x = (xd .- x_tangential) ./ radial
        y = (yd .- y_tangential) ./ radial
    end

    R_matrix = Float64.(R)
    R_rel_matrix = Float64.(R_rel)
    translation = Float64.(t)
    relative_translation = Float64.(t_rel)
    camera_ray = vcat(reshape(x, 1, :), reshape(y, 1, :), ones(1, length(x)))
    camera_center_marker = -(R_matrix' * translation)
    ray_marker = R_matrix' * camera_ray
    model_plane_point = relative_translation .+ R_rel_matrix[:, 3] .* Float64(plane_coordinate)
    model_plane_normal = R_rel_matrix[:, 3]
    ray_scale = dot(model_plane_normal, model_plane_point .- camera_center_marker) ./
                vec(model_plane_normal' * ray_marker)
    X_marker = camera_center_marker .+ ray_marker .* reshape(ray_scale, 1, :)
    X_model = R_rel_matrix' * (X_marker .- relative_translation)

    return X_model[1, :], X_model[2, :], X_model[3, :]
end

function pixel_to_model(
    u::Real,
    v::Real,
    R::AbstractMatrix,
    t::AbstractVector,
    K::AbstractMatrix,
    dist_coeffs::AbstractVector,
    R_rel::AbstractMatrix,
    t_rel::AbstractVector;
    kwargs...,
)
    r, z, s = pixel_to_model(
        [u], [v], R, t, K, dist_coeffs, R_rel, t_rel;
        kwargs...,
    )
    return r[1], z[1], s[1]
end




function getRrel_fn(ga::Float64, gb::Float64, gw::Float64)
    """
    Get the relative rotation matrix from Aruco marker to droplet

    arg:
    - ga    : rotation angle around i_w [degrees]
    - gb    : rotation angle around j_w [degrees]
    - gw    : rotation angle around k_w [degrees]

    return:
    - R_rel : relative rotation matrix
    - t_rel : relative translation vector
    """
        
    Ri = Float64.([
        1      0       0;
        0  cosd(ga) -sind(ga);
        0  sind(ga)  cosd(ga)
    ])

    Rj = Float64.([
        cosd(gb) 0.0 sind(gb);
        0.0     1.0 0.0;
        -sind(gb) 0.0 cosd(gb)
    ])

    Rk = Float64.([
        cosd(gw) -sind(gw) 0.0;
        sind(gw)  cosd(gw) 0.0;
        0.0      0.0     1.0

    ])

    return Rk * Rj * Ri
end

function R_to_euler_zyx(R)

    pitch = asind(-R[3,1])

    roll = atand(R[3,2], R[3,3])

    yaw = atand(R[2,1], R[1,1])

    return roll, pitch, yaw
end



function save_calibration_results(calib_results, output_path::String)
    if calib_results === nothing
        println("No calibration results available to save.")
        return false
    end

    cam_matrix  = calib_results["camera_matrix"]
    dist_coeffs = calib_results["dist_coeffs"]
    rms         = calib_results["rms"]

    @save output_path cam_matrix dist_coeffs rms
    println("Saved calibration JLD2 file at: $output_path")
    println("RMS reprojection error: $rms")
    return true
end


function plot_intrinsic_distributions(result; show_plots::Bool=true, save_fig::Bool=false, out_dir::String=".")
    """
    Plot distributions of camera intrinsic matrix components and distortion coefficients
    from the `result` produced by `uncertainty_camera_intrinsics`.

    `result` is expected to be a Dict with keys `k_samples` (array of 3x3 lists)
    and `dist_samples` (array of 5-element lists).
    """

    k_samples = result["k_samples"]
    dist_samples = result["dist_samples"]

    if isempty(k_samples) && isempty(dist_samples)
        @warn "No samples to plot"
        return nothing
    end

    # Flatten k_samples into arrays for the six unique entries we care about
    # fx = K[1,1], fy = K[2,2], cx = K[1,3], cy = K[2,3]
    fx = [k[0][0] for k in k_samples]
    fy = [k[1][1] for k in k_samples]
    cx = [k[0][2] for k in k_samples]
    cy = [k[1][2] for k in k_samples]

    # dist: k1, k2, p1, p2, k3
    k1 = [d[0] for d in dist_samples]
    k2 = [d[1] for d in dist_samples]
    p1 = [d[2] for d in dist_samples]
    p2 = [d[3] for d in dist_samples]
    k3 = [d[4] for d in dist_samples]

    plt = plot(layout = (3,2), size=(900,900))
    histogram!(plt[1], fx, bins=30, title="fx", legend=false)
    histogram!(plt[2], fy, bins=30, title="fy", legend=false)
    histogram!(plt[3], cx, bins=30, title="cx", legend=false)
    histogram!(plt[4], cy, bins=30, title="cy", legend=false)
    histogram!(plt[5], k1, bins=30, title="k1", legend=false)
    histogram!(plt[6], k2, bins=30, title="k2", legend=false)

    if save_fig
        ensure_directory(out_dir)
        png_path = joinpath(out_dir, "intrinsic_distributions.png")
        savefig(plt, png_path)
    end

    if show_plots
        display(plt)
    end

    return plt
end

module IntrinsicUncertainty
using ..Camera: ensure_directory
using FileIO, JLD2, Plots, PyCall, Random, Statistics

include(joinpath(@__DIR__, "camera_uncertainties.jl"))
end

end