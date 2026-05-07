using PyCall
using FileIO
using Images

function ensure_directory(path::String)
    if !isdir(path)
        mkpath(path)
    end
end

function estimate_homography_single(img_path::String, K_path::String, dist_path::String, 
                                   marker_length_mm::Float64=5.0;
                                   px_per_mm::Int=50, bounds::Int=10, 
                                   output_path::String="homography_output.jpg")
    """
    Estimate homography for a single image with detected ArUco marker.
    Creates an orthographic projection of the marker plane.
    
    Args:
        img_path: Path to input image with detected marker
        K_path: Path to camera matrix K (numpy .npy file)
        dist_path: Path to distortion coefficients (numpy .npy file)
        marker_length_mm: Physical size of marker in mm
        px_per_mm: Output resolution (pixels per mm)
        bounds: Canvas expansion factor (bounds=1 -> only marker, bounds=10 -> 10x marker size)
        output_path: Path to save homography-warped image
    
    Returns:
        Dict with homography results
    """
    
    result = py"""
import cv2, numpy as np, sys

# --- Load parameters ---
K = np.load($K_path)
dist = np.load($dist_path)

marker_length_mm = $marker_length_mm
px_per_mm = $px_per_mm
bounds = $bounds

# --- Read experimental image ---
img = cv2.imread($img_path)

if img is None:
    print(f"Could not read {$img_path}")
    sys.stdout.flush()
    None
else:
    # --- Ensure 3-channel BGR ---
    if img.ndim == 2:
        img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
    elif img.shape[2] == 4:
        img = cv2.cvtColor(img, cv2.COLOR_RGBA2BGR)

    # --- Convert to gray ---
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

    # --- Define ArUco dictionary + parameters ---
    aruco_dict = cv2.aruco.getPredefinedDictionary(cv2.aruco.DICT_4X4_50)
    params = cv2.aruco.DetectorParameters()

    # --- Create ArUco detector ---
    detector = cv2.aruco.ArucoDetector(aruco_dict, params)

    # --- Detect ArUco marker ---
    corners, ids, rejected = detector.detectMarkers(gray)
    print("Detected marker IDs:", ids)
    sys.stdout.flush()

    # --- If a marker was detected in the image ---
    if ids is not None:
        # --- Extract detected marker corners & reorganize ---
        # corners[0] has shape (1, 4, 2); reshape to (4,2)
        # OpenCV ArUco returns points in this order:
        #   0: top-left, 1: top-right, 2: bottom-right, 3: bottom-left
        corners = corners[0].reshape((4, 2))

        # --- Define points in world coordinates [mm] ---
        # (order needs to match order of corners)
        world_pts = np.array([
            [0.0, 0.0],
            [marker_length_mm, 0.0],
            [marker_length_mm, marker_length_mm],
            [0.0, marker_length_mm]
        ], dtype=np.float32)

        img_pts = corners.astype(np.float32)

        # --- Compute output image dimensions ---
        out_w = int(marker_length_mm * px_per_mm)
        out_h = int(marker_length_mm * px_per_mm)

        # --- Define points in destination image coordinates [px] ---
        # These are the positions where the marker corners should land after warping
        dst_pts = np.array([
            [0, 0],
            [out_w - 1, 0],
            [out_w - 1, out_h - 1],
            [0, out_h - 1]
        ], dtype=np.float32)

        # --- Compute homography mapping image points -> orthographic marker plane ---
        H, status = cv2.findHomography(img_pts, dst_pts)

        # --- Define overall output image size ---
        canvas_w = int(marker_length_mm * px_per_mm * bounds)
        canvas_h = int(marker_length_mm * px_per_mm * bounds)

        # --- Place the marker in the centre of the output image ---
        tx = (canvas_w - out_w) / 2.0
        ty = (canvas_h - out_h) / 2.0
        T = np.array([[1, 0, tx], [0, 1, ty], [0, 0, 1]], dtype=np.float64)
        H_centered = T.dot(H)

        # --- Warp image into orthographic projection based on homography ---
        img_hom = cv2.warpPerspective(img, H_centered, (canvas_w, canvas_h), flags=cv2.INTER_LINEAR)

        # --- Sanity check: measure the pixel size of the warped marker ---
        mapped = cv2.perspectiveTransform(img_pts.reshape(-1, 1, 2), H_centered).reshape(-1, 2)
        p0, p1 = mapped[0], mapped[1]
        pixel_length = np.linalg.norm(p1 - p0)
        px_per_mm_measured = pixel_length / marker_length_mm

        print(f"Px/mm set: {px_per_mm}; Measured px/mm in result: {px_per_mm_measured:.2f}")
        print(f"Canvas size: {canvas_w}x{canvas_h} px")
        print(f"Homography matrix:\n{H_centered}")
        sys.stdout.flush()

        # --- Save warped image ---
        cv2.imwrite($output_path, img_hom)
        print(f"Homography image saved to {$output_path}")
        sys.stdout.flush()

        {
            "success": True,
            "marker_ids": ids.flatten().tolist(),
            "canvas_size": (canvas_w, canvas_h),
            "marker_size_px": (out_w, out_h),
            "px_per_mm_measured": float(px_per_mm_measured),
            "H_matrix": H_centered.tolist(),
            "output_path": $output_path
        }
    else:
        print("No markers detected")
        sys.stdout.flush()
        {"success": False, "message": "No markers detected"}
"""

    return result
end

function estimate_homography_batch(img_dir::String, out_dir::String, K_path::String, dist_path::String,
                                  marker_length_mm::Float64=5.0;
                                  px_per_mm::Int=50, bounds::Int=10)
    """
    Estimate homography for all images with detected ArUco markers in a directory.
    
    Args:
        img_dir: Directory containing images with detected markers
        out_dir: Directory to save homography-warped images
        K_path: Path to camera matrix K (numpy .npy file)
        dist_path: Path to distortion coefficients (numpy .npy file)
        marker_length_mm: Physical size of marker in mm
        px_per_mm: Output resolution (pixels per mm)
        bounds: Canvas expansion factor
    
    Returns:
        Dict with batch processing results
    """
    
    ensure_directory(out_dir)
    
    result = py"""
import cv2, numpy as np, os, glob, sys

# --- Load parameters ---
K = np.load($K_path)
dist = np.load($dist_path)

marker_length_mm = $marker_length_mm
px_per_mm = $px_per_mm
bounds = $bounds

# --- Prep output directory ---
os.makedirs($out_dir, exist_ok=True)

# --- Read experimental images ---
img_files = glob.glob(os.path.join($img_dir, "*.JPG")) + glob.glob(os.path.join($img_dir, "*.jpg"))

homography_results = {
    "processed_count": 0,
    "error_count": 0,
    "files": []
}

for img_path in img_files:
    img = cv2.imread(img_path)

    if img is None:
        print(f"Could not read {img_path}")
        sys.stdout.flush()
        homography_results["error_count"] += 1
        continue

    # --- Ensure 3-channel BGR ---
    if img.ndim == 2:
        img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
    elif img.shape[2] == 4:
        img = cv2.cvtColor(img, cv2.COLOR_RGBA2BGR)

    # --- Convert to gray ---
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

    # --- Define ArUco dictionary + parameters ---
    aruco_dict = cv2.aruco.getPredefinedDictionary(cv2.aruco.DICT_4X4_50)
    params = cv2.aruco.DetectorParameters()

    # --- Create ArUco detector ---
    detector = cv2.aruco.ArucoDetector(aruco_dict, params)

    # --- Detect ArUco marker ---
    corners, ids, rejected = detector.detectMarkers(gray)

    # --- If marker detected ---
    if ids is not None:
        corners = corners[0].reshape((4, 2))

        # --- Define points in world coordinates [mm] ---
        world_pts = np.array([
            [0.0, 0.0],
            [marker_length_mm, 0.0],
            [marker_length_mm, marker_length_mm],
            [0.0, marker_length_mm]
        ], dtype=np.float32)

        img_pts = corners.astype(np.float32)

        out_w = int(marker_length_mm * px_per_mm)
        out_h = int(marker_length_mm * px_per_mm)

        dst_pts = np.array([
            [0, 0],
            [out_w - 1, 0],
            [out_w - 1, out_h - 1],
            [0, out_h - 1]
        ], dtype=np.float32)

        H, status = cv2.findHomography(img_pts, dst_pts)

        canvas_w = int(marker_length_mm * px_per_mm * bounds)
        canvas_h = int(marker_length_mm * px_per_mm * bounds)

        tx = (canvas_w - out_w) / 2.0
        ty = (canvas_h - out_h) / 2.0
        T = np.array([[1, 0, tx], [0, 1, ty], [0, 0, 1]], dtype=np.float64)
        H_centered = T.dot(H)

        img_hom = cv2.warpPerspective(img, H_centered, (canvas_w, canvas_h), flags=cv2.INTER_LINEAR)

        out_img_path = os.path.join($out_dir, os.path.basename(img_path))
        cv2.imwrite(out_img_path, img_hom)

        homography_results["processed_count"] += 1
        homography_results["files"].append({
            "filename": os.path.basename(img_path),
            "marker_ids": ids.flatten().tolist(),
            "canvas_size": (canvas_w, canvas_h)
        })

        print(f"Homography computed for {os.path.basename(img_path)}")
        sys.stdout.flush()

    else:
        print(f"No markers detected in {os.path.basename(img_path)}")
        sys.stdout.flush()
        homography_results["error_count"] += 1

print(f"\nHomography images saved to {$out_dir}")
sys.stdout.flush()

homography_results
"""

    return result
end

function print_usage()
    println("Usage: julia EstimateHomographyFromNotebook.jl [command] [args...]")
    println()
    println("Commands:")
    println("  single [img_path] [K_path] [dist_path] [output_path] [marker_length_mm] [px_per_mm] [bounds]")
    println("    Estimate homography for a single image")
    println()
    println("  batch [img_dir] [out_dir] [K_path] [dist_path] [marker_length_mm] [px_per_mm] [bounds]")
    println("    Estimate homography for all images in a directory")
    println()
    println("Parameters:")
    println("  marker_length_mm : Physical marker size (default: 5.0)")
    println("  px_per_mm        : Output resolution in pixels per mm (default: 50)")
    println("  bounds           : Canvas expansion factor (default: 10)")
    println()
    println("Examples:")
    println("  julia EstimateHomographyFromNotebook.jl single \"image.jpg\" K.npy dist.npy \"output.jpg\"")
    println("  julia EstimateHomographyFromNotebook.jl batch \"./detected\" \"./homography\" K.npy dist.npy 5.0 50 10")
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) == 0 || ARGS[1] in ["-h", "--help"]
        print_usage()
        exit()
    end

    command = ARGS[1]

    if command == "single" && length(ARGS) >= 4
        img_path = ARGS[2]
        K_path = ARGS[3]
        dist_path = ARGS[4]
        output_path = length(ARGS) >= 5 ? ARGS[5] : "homography_output.jpg"
        marker_length_mm = length(ARGS) >= 6 ? parse(Float64, ARGS[6]) : 5.0
        px_per_mm = length(ARGS) >= 7 ? parse(Int, ARGS[7]) : 50
        bounds = length(ARGS) >= 8 ? parse(Int, ARGS[8]) : 10

        println("Homography Estimation - Single Image Mode")
        println("=" * 50)
        println("Input image: $img_path")
        println("Output image: $output_path")
        println("Camera matrix: $K_path")
        println("Distortion coeffs: $dist_path")
        println("Marker length: $marker_length_mm mm")
        println("Px per mm: $px_per_mm")
        println("Canvas bounds: $bounds")
        println()

        if !isfile(img_path)
            println("✗ Error: Image file not found: $img_path")
            exit(1)
        end
        if !isfile(K_path)
            println("✗ Error: Camera matrix file not found: $K_path")
            exit(1)
        end
        if !isfile(dist_path)
            println("✗ Error: Distortion coefficients file not found: $dist_path")
            exit(1)
        end

        result = estimate_homography_single(img_path, K_path, dist_path, marker_length_mm;
                                           px_per_mm=px_per_mm, bounds=bounds,
                                           output_path=output_path)

        if result !== nothing && result["success"] == true
            println("\n=== Homography Results ===")
            println("✓ Marker IDs: $(result["marker_ids"])")
            println("Canvas size: $(result["canvas_size"]) px")
            println("Marker size: $(result["marker_size_px"]) px")
            println("Measured px/mm: $(round(result["px_per_mm_measured"]; digits=2))")
        else
            println("\n✗ Failed to estimate homography")
            exit(1)
        end

    elseif command == "batch" && length(ARGS) >= 5
        img_dir = ARGS[2]
        out_dir = ARGS[3]
        K_path = ARGS[4]
        dist_path = ARGS[5]
        marker_length_mm = length(ARGS) >= 6 ? parse(Float64, ARGS[6]) : 5.0
        px_per_mm = length(ARGS) >= 7 ? parse(Int, ARGS[7]) : 50
        bounds = length(ARGS) >= 8 ? parse(Int, ARGS[8]) : 10

        println("Homography Estimation - Batch Mode")
        println("=" * 50)
        println("Input directory: $img_dir")
        println("Output directory: $out_dir")
        println("Camera matrix: $K_path")
        println("Distortion coeffs: $dist_path")
        println("Marker length: $marker_length_mm mm")
        println("Px per mm: $px_per_mm")
        println("Canvas bounds: $bounds")
        println()

        if !isdir(img_dir)
            println("✗ Error: Image directory not found: $img_dir")
            exit(1)
        end
        if !isfile(K_path)
            println("✗ Error: Camera matrix file not found: $K_path")
            exit(1)
        end
        if !isfile(dist_path)
            println("✗ Error: Distortion coefficients file not found: $dist_path")
            exit(1)
        end

        results = estimate_homography_batch(img_dir, out_dir, K_path, dist_path, marker_length_mm;
                                           px_per_mm=px_per_mm, bounds=bounds)

        if results !== nothing
            println("\n=== Batch Processing Summary ===")
            println("✓ Successfully processed: $(results["processed_count"]) images")
            if results["error_count"] > 0
                println("✗ Errors: $(results["error_count"]) images")
            end
        end

    else
        println("Unknown command or missing arguments: $command")
        println()
        print_usage()
        exit(1)
    end
end
