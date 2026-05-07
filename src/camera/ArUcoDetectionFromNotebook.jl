using PyCall
using FileIO
using Images

function ensure_directory(path::String)
    if !isdir(path)
        mkpath(path)
    end
end

function detect_aruco_markers(img_dir::String, out_dir::String, K_path::String, dist_path::String, marker_length_mm::Float64=5.0)
    ensure_directory(out_dir)

    py"""
import cv2, numpy as np, os, glob, sys

# --- Load parameters ---
K = np.load($K_path)
dist = np.load($dist_path)
marker_length_mm = $marker_length_mm

# --- Prep output directory ---
os.makedirs($out_dir, exist_ok=True)

# --- Read experimental images ---
img_files = glob.glob(os.path.join($img_dir, "*.JPG")) + glob.glob(os.path.join($img_dir, "*.jpg"))

detection_results = []

for img_path in img_files:
    img = cv2.imread(img_path)

    if img is None:
        print(f"Could not read {img_path}")
        sys.stdout.flush()
        continue

    # --- Ensure 3-channel BGR ---
    if img.ndim == 2:
        img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
    elif img.shape[2] == 4:
        img = cv2.cvtColor(img, cv2.COLOR_RGBA2BGR)

    # --- Convert to gray ---
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

    # --- Define ArUco dictionary + parameters used ---
    # Make sure you use the same dictionary that was used to generate the marker
    aruco_dict = cv2.aruco.getPredefinedDictionary(cv2.aruco.DICT_4X4_50)
    params = cv2.aruco.DetectorParameters()

    # --- Create ArUco detector ---
    detector = cv2.aruco.ArucoDetector(aruco_dict, params)

    # --- Detect ArUco marker ---
    corners, ids, rejected = detector.detectMarkers(gray)
    print("Detected marker IDs:", ids)

    # --- If a marker was detected in the image ---
    if ids is not None:
        # --- Estimate pose ---
        rvecs, tvecs, obj_points = cv2.aruco.estimatePoseSingleMarkers(corners, marker_length_mm, K, dist)
        rvec, tvec = rvecs[0], tvecs[0]
        marker_dist = tvec[0][2]

        # --- Pixel/mm scaling ---
        img_points, _ = cv2.projectPoints(obj_points, rvec, tvec, K, dist)
        img_points = img_points.reshape(-1, 2)
        
        # Distance between two corners of the marker
        px_dist = np.linalg.norm(img_points[0] - img_points[1])
        mm_dist = marker_length_mm
        
        # Configure px/mm scaling
        px_per_mm = px_dist / mm_dist
        mm_per_px = mm_dist / px_dist

        print(f"Scale: {px_per_mm: .2f} px/mm ({mm_per_px:.3f} mm/px)")
        print(f"Marker distance from camera: {marker_dist:.2f} mm\n")
        sys.stdout.flush()

        # --- Save images where a marker was detected ---
        out_img_path = os.path.join($out_dir, os.path.basename(img_path))
        cv2.imwrite(out_img_path, img)

        # Store result
        detection_results.append({
            "filename": os.path.basename(img_path),
            "marker_ids": ids.flatten().tolist() if ids is not None else [],
            "marker_dist_mm": float(marker_dist),
            "px_per_mm": float(px_per_mm),
            "mm_per_px": float(mm_per_px),
        })

        print("Image saved to", $out_dir)
        sys.stdout.flush()

    # Else, no marker was detected
    else:
        print("No markers detected in", os.path.basename(img_path))
        sys.stdout.flush()

print("\nImages with detected markers saved to", $out_dir)
detection_results
"""

    return py"detection_results"
end

function process_single_image(img_path::String, K_path::String, dist_path::String, marker_length_mm::Float64=5.0)
    """Process a single image and display results"""
    
    result = py"""
import cv2, numpy as np, sys

# --- Load parameters ---
K = np.load($K_path)
dist = np.load($dist_path)
marker_length_mm = $marker_length_mm

# --- Read experimental image ---
img = cv2.imread($img_path)

if img is None:
    print("Image not loaded. Please check the image path.")
    sys.stdout.flush()
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
    print("Detected marker IDs:", ids, "\n")
    sys.stdout.flush()

    if ids is not None:
        # Draw detected marker
        marker = cv2.aruco.drawDetectedMarkers(img, corners, ids)
        cv2.imwrite("marker_detected.jpg", marker)

        # --- Estimate pose ---
        rvecs, tvecs, obj_points = cv2.aruco.estimatePoseSingleMarkers(corners, marker_length_mm, K, dist)
        rvec, tvec = rvecs[0], tvecs[0]

        # --- Pixel/mm scaling ---
        img_points, _ = cv2.projectPoints(obj_points, rvec, tvec, K, dist)
        img_points = img_points.reshape(-1, 2)

        px_dist = np.linalg.norm(img_points[0] - img_points[1])
        mm_dist = marker_length_mm
        px_per_mm = px_dist / mm_dist
        mm_per_px = mm_dist / px_dist

        print(f"Scale: {px_per_mm: .2f} px/mm ({mm_per_px:.3f} mm/px)")
        print(f"Marker distance from camera: {tvec[0][2]:.2f} mm\n")

        # --- Focal length ---
        W_sensor = 6.17  # width of physical camera sensor [mm]
        H_sensor = 4.58  # height of physical camera sensor [mm]
        N_x = img.shape[1]  # width of image [pixels]
        N_y = img.shape[0]  # height of image [pixels]

        f_x = K[0][0]  # focal length in x [px]
        f_y = K[1][1]  # focal length in y [px]

        f_mm_x = f_x / (N_x / W_sensor)  # focal length in x [mm]
        f_mm_y = f_y / (N_y / H_sensor)  # focal length in y [mm]

        f_camera = 10  # focal length of camera in image metadata

        err_x = abs(f_camera - f_mm_x)
        err_y = abs(f_camera - f_mm_y)

        print("FOCAL LENGTH INFORMATION")
        print(f"PX: fx = {f_x:.2f} px, fy = {f_y:.2f} px")
        print(f"MM: fx = {f_mm_x:.2f} mm, fy = {f_mm_y:.2f} mm")
        print(f"Error: |f-fx| = {err_x:.2f} mm, |f-fy| = {err_y:.2f} mm")
        sys.stdout.flush()

        {
            "success": True,
            "marker_ids": ids.flatten().tolist(),
            "marker_dist_mm": float(tvec[0][2]),
            "px_per_mm": float(px_per_mm),
            "mm_per_px": float(mm_per_px),
            "f_x": float(f_x),
            "f_y": float(f_y),
            "f_mm_x": float(f_mm_x),
            "f_mm_y": float(f_mm_y),
        }
    else:
        print("No markers detected")
        sys.stdout.flush()
        {"success": False}
"""

    return result
end

function print_usage()
    println("Usage: julia ArUcoDetectionFromNotebook.jl [command] [args...]")
    println()
    println("Commands:")
    println("  batch [img_dir] [out_dir] [K_path] [dist_path] [marker_length_mm]")
    println("    Detect ArUco markers in all images in a directory")
    println()
    println("  single [img_path] [K_path] [dist_path] [marker_length_mm]")
    println("    Process a single image with focal length information")
    println()
    println("Example:")
    println("  julia ArUcoDetectionFromNotebook.jl batch \"./images\" \"./detected\" \"K.npy\" \"dist.npy\" 5.0")
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) == 0 || ARGS[1] in ["-h", "--help"]
        print_usage()
        exit()
    end

    command = ARGS[1]

    if command == "batch" && length(ARGS) >= 5
        img_dir = ARGS[2]
        out_dir = ARGS[3]
        K_path = ARGS[4]
        dist_path = ARGS[5]
        marker_length_mm = length(ARGS) >= 6 ? parse(Float64, ARGS[6]) : 5.0

        println("ArUco marker detection - Batch mode")
        println("  image directory: $img_dir")
        println("  output directory: $out_dir")
        println("  K path: $K_path")
        println("  dist path: $dist_path")
        println("  marker length: $marker_length_mm mm")
        println()

        results = detect_aruco_markers(img_dir, out_dir, K_path, dist_path, marker_length_mm)

        if results !== nothing && length(results) > 0
            println("\n=== Detection Summary ===")
            for result in results
                println("File: $(result["filename"])")
                println("  Marker IDs: $(result["marker_ids"])")
                println("  Distance: $(result["marker_dist_mm"]:.2f) mm")
                println("  Scale: $(result["px_per_mm"]:.2f) px/mm")
            end
        else
            println("No markers detected in any images.")
        end

    elseif command == "single" && length(ARGS) >= 4
        img_path = ARGS[2]
        K_path = ARGS[3]
        dist_path = ARGS[4]
        marker_length_mm = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 5.0

        println("ArUco marker detection - Single image mode")
        println("  image: $img_path")
        println("  K path: $K_path")
        println("  dist path: $dist_path")
        println("  marker length: $marker_length_mm mm")
        println()

        result = process_single_image(img_path, K_path, dist_path, marker_length_mm)

        if result !== nothing && result["success"] == true
            println("\n=== Detection Results ===")
            println("Marker IDs: $(result["marker_ids"])")
            println("Distance: $(result["marker_dist_mm"]:.2f} mm")
            println("Scale: $(result["px_per_mm"]:.2f} px/mm")
        else
            println("No markers detected.")
        end

    else
        println("Unknown command or missing arguments: $command")
        println()
        print_usage()
        exit(1)
    end
end
