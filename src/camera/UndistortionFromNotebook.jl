using PyCall

function ensure_directory(path::String)
    if !isdir(path)
        mkpath(path)
    end
end

function undistort_images(img_dir::String, out_dir::String, K_path::String, dist_path::String; 
                         draw_roi::Bool=false, verbose::Bool=true)
    """
    Undistort a batch of experimental images using calibration parameters.
    
    Args:
        img_dir: Directory containing experimental images
        out_dir: Directory to save undistorted images
        K_path: Path to camera matrix K (numpy .npy file)
        dist_path: Path to distortion coefficients (numpy .npy file)
        draw_roi: Whether to draw region of interest rectangle on output images
        verbose: Whether to print status messages
    
    Returns:
        Dict with processing results
    """
    
    ensure_directory(out_dir)
    
    result = py"""
import cv2, numpy as np, os, glob, sys

# --- Load parameters ---
K = np.load($K_path)
dist = np.load($dist_path)

# --- Prep output directory ---
os.makedirs($out_dir, exist_ok=True)

# --- Read experimental images ---
img_files = glob.glob(os.path.join($img_dir, "*.JPG")) + glob.glob(os.path.join($img_dir, "*.jpg"))

undistort_results = {
    "processed_count": 0,
    "error_count": 0,
    "files": []
}

for img_path in img_files:
    img = cv2.imread(img_path)
    if img is None:
        print(f"Could not read {img_path}")
        sys.stdout.flush()
        undistort_results["error_count"] += 1
        continue
    
    h, w = img.shape[:2]
    size = (w, h)

    # --- Compute optimal new camera matrix ---
    K_new, roi = cv2.getOptimalNewCameraMatrix(K, dist, size, 1.0, size)
    
    # --- Undistort image ---
    ud_img = cv2.undistort(img, K, dist, None, K_new)
    
    # --- Optional: draw ROI for visualization ---
    if $draw_roi:
        x, y, w_roi, h_roi = roi
        cv2.rectangle(ud_img, (x,y), (x+w_roi, y+h_roi), (0, 255, 0), 2)
    
    # --- Save ---
    out_img_path = os.path.join($out_dir, os.path.basename(img_path))
    cv2.imwrite(out_img_path, ud_img)
    
    undistort_results["processed_count"] += 1
    undistort_results["files"].append({
        "filename": os.path.basename(img_path),
        "K_new": K_new.tolist(),
        "roi": list(roi)
    })
    
    if $verbose:
        print(f"Undistorted: {os.path.basename(img_path)}")
        sys.stdout.flush()

if $verbose:
    print(f"\nUndistorted images saved to {$out_dir}")
    sys.stdout.flush()

undistort_results
"""
    
    return result
end

function undistort_single_image(img_path::String, K_path::String, dist_path::String; 
                               output_path::String="undistorted_output.jpg",
                               draw_roi::Bool=false)
    """
    Undistort a single experimental image.
    
    Args:
        img_path: Path to input image
        K_path: Path to camera matrix K (numpy .npy file)
        dist_path: Path to distortion coefficients (numpy .npy file)
        output_path: Path to save undistorted image
        draw_roi: Whether to draw region of interest rectangle
    
    Returns:
        Dict with undistortion results
    """
    
    result = py"""
import cv2, numpy as np, sys

# --- Load parameters ---
K = np.load($K_path)
dist = np.load($dist_path)

# --- Read experimental image ---
img = cv2.imread($img_path)

if img is None:
    print(f"Could not read {$img_path}")
    sys.stdout.flush()
    None
else:
    h, w = img.shape[:2]
    size = (w, h)

    # --- Compute optimal new camera matrix ---
    K_new, roi = cv2.getOptimalNewCameraMatrix(K, dist, size, 1.0, size)
    
    # --- Undistort image ---
    ud_img = cv2.undistort(img, K, dist, None, K_new)
    
    # --- Optional: draw ROI for visualization ---
    if $draw_roi:
        x, y, w_roi, h_roi = roi
        cv2.rectangle(ud_img, (x,y), (x+w_roi, y+h_roi), (0, 255, 0), 2)
    
    # --- Save ---
    cv2.imwrite($output_path, ud_img)
    
    print(f"Undistorted image saved to {$output_path}")
    print(f"Image size: {w}x{h}")
    print(f"ROI: {roi}")
    sys.stdout.flush()
    
    {
        "success": True,
        "image_size": (w, h),
        "K_new": K_new.tolist(),
        "roi": list(roi),
        "output_path": $output_path
    }
"""
    
    return result
end

function print_usage()
    println("Usage: julia UndistortionFromNotebook.jl [command] [args...]")
    println()
    println("Commands:")
    println("  batch [img_dir] [out_dir] [K_path] [dist_path]")
    println("    Undistort all images in a directory")
    println()
    println("  single [img_path] [K_path] [dist_path] [output_path]")
    println("    Undistort a single image")
    println()
    println("Example:")
    println("  julia UndistortionFromNotebook.jl batch \"./images\" \"./undistorted\" K.npy dist.npy")
    println("  julia UndistortionFromNotebook.jl single \"image.jpg\" K.npy dist.npy \"undistorted.jpg\"")
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
        draw_roi = length(ARGS) >= 6 && ARGS[6] in ["true", "True", "1"]

        println("Image Undistortion - Batch Mode")
        println("=" * 50)
        println("Input directory: $img_dir")
        println("Output directory: $out_dir")
        println("Camera matrix: $K_path")
        println("Distortion coeffs: $dist_path")
        println("Draw ROI: $draw_roi")
        println()

        if !isfile(K_path)
            println("✗ Error: Camera matrix file not found: $K_path")
            exit(1)
        end
        
        if !isfile(dist_path)
            println("✗ Error: Distortion coefficients file not found: $dist_path")
            exit(1)
        end
        
        if !isdir(img_dir)
            println("✗ Error: Image directory not found: $img_dir")
            exit(1)
        end

        results = undistort_images(img_dir, out_dir, K_path, dist_path; 
                                  draw_roi=draw_roi, verbose=true)

        if results !== nothing
            println("\n=== Processing Summary ===")
            println("✓ Successfully processed: $(results["processed_count"]) images")
            if results["error_count"] > 0
                println("✗ Errors: $(results["error_count"]) images")
            end
        end

    elseif command == "single" && length(ARGS) >= 4
        img_path = ARGS[2]
        K_path = ARGS[3]
        dist_path = ARGS[4]
        output_path = length(ARGS) >= 5 ? ARGS[5] : "undistorted_output.jpg"
        draw_roi = length(ARGS) >= 6 && ARGS[6] in ["true", "True", "1"]

        println("Image Undistortion - Single Image Mode")
        println("=" * 50)
        println("Input image: $img_path")
        println("Output image: $output_path")
        println("Camera matrix: $K_path")
        println("Distortion coeffs: $dist_path")
        println("Draw ROI: $draw_roi")
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

        result = undistort_single_image(img_path, K_path, dist_path; 
                                       output_path=output_path, 
                                       draw_roi=draw_roi)

        if result !== nothing && result["success"] == true
            println("\n✓ Image successfully undistorted!")
        else
            println("\n✗ Failed to undistort image")
            exit(1)
        end

    else
        println("Unknown command or missing arguments: $command")
        println()
        print_usage()
        exit(1)
    end
end
