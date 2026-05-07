using PyCall
using FileIO
using JLD2

function ensure_directory(path::String)
    if !isdir(path)
        mkpath(path)
    end
end

function calibrate_camera(img_dir::String, out_dir::String; checkerboard_size::Tuple{Int,Int}=(13, 9), square_size::Float64=0.021)
    """
    Calibrate camera from checkerboard images.
    
    Detects checkerboard corners in calibration images, estimates camera matrix K 
    and distortion coefficients, and saves results as .jld2 files.
    
    Args:
        img_dir:            Directory containing checkerboard calibration images
        out_dir:            Directory to save detected images and calibration files
        checkerboard_size:  Tuple of (horizontal_corners, vertical_corners)
        square_size:        Physical size of checkerboard squares in meters
    
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
img_files = glob.glob(os.path.join($img_dir, "*.JPG")) + glob.glob(os.path.join($img_dir, "*.jpg"))

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
    ret, corners    = cv2.findChessboardCorners(gray, checker, flags=flags)

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

# --- Calibrate camera if corners were found ---
if len(objpoints) > 0:
    h, w                            = gray.shape[:2]
    ret, mtx, dist, rvecs, tvecs    = cv2.calibrateCamera(objpoints, imgpoints, (w, h), None, None)

    # --- Save calibration matrices ---
    np.save(os.path.join($out_dir, "camera_matrix.npy"), mtx)
    np.save(os.path.join($out_dir, "dist_coeffs.npy"), dist)

    {
        "rms": float(ret),
        "camera_matrix": mtx.tolist(),
        "dist_coeffs": dist.tolist(),
        "detected_paths": list(saved_paths),
    }
else:
    print("No corners detected; calibration skipped")
    None
"""

    return result
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

function print_usage()
    println("Usage: julia CameraCalibrationFromNotebook.jl [img_dir] [out_dir] [set_id]")
    println("  img_dir   : directory containing calibration images")
    println("  out_dir   : directory where detected images and results will be saved")
    println("  set_id    : optional integer suffix for the output JLD2 file")
    println("If no arguments are provided, default paths in the script are used.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    # --- Configuration Defaults ---
    dir_img_calibration         = "data/calibration_camera"
    default_img_dir             = joinpath(dir_img_calibration, "images")
    default_out_dir             = joinpath(dir_img_calibration, "detected")
    default_checkerboard_size   = (13, 9)
    default_square_size         = 0.021
    default_set_id              = 1

    if length(ARGS) == 1 && ARGS[1] in ["-h", "--help"]
        print_usage()
        exit()
    end

    img_dir = length(ARGS) >= 1 ? ARGS[1] : default_img_dir
    out_dir = length(ARGS) >= 2 ? ARGS[2] : default_out_dir
    set_id  = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : default_set_id

    println("Camera calibration starting")
    println("  image directory: $img_dir")
    println("  output directory: $out_dir")
    println("  checkerboard size: $default_checkerboard_size")
    println("  square size: $default_square_size")
    println()

    calib_results = calibrate_camera(img_dir, out_dir, default_checkerboard_size, default_square_size)

    if calib_results !== nothing
        jld2_path = joinpath(out_dir, "calibration_set$(set_id).jld2")
        save_calibration_results(calib_results, jld2_path)
    else
        println("Calibration did not produce results.")
    end
end
