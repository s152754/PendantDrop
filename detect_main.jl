using JLD2
using PyCall
using FileIO
using Images
using Plots

function Temp_detect(img_file)
    py"""
    import cv2, numpy as np, os, sys

    img_path            = $img_file

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
    detector    = cv2.aruco.ArucoDetector(aruco_dict, params)

    corners, ids, rejected = detector.detectMarkers(gray)
    if ids is None:
        raise ValueError(f"No markers detected in {img_path}")

    """
    return Array(py"ids")
end

# Folder containing the images to batch process
img_folder = "test_images"  # Change this to your image folder path
detected_folder = joinpath(img_folder, "detected")
not_detected_folder = joinpath(img_folder, "undetected")

# Create output folders if they don't exist
mkpath(detected_folder)
mkpath(not_detected_folder)

# Get all JPG files in the folder
img_files = filter(f -> lowercase(f)[end-3:end] == ".jpg", readdir(img_folder))

println("Processing $(length(img_files)) JPG files from: $img_folder")
println()

detected = String[]
not_detected = String[]

for img_file in img_files
    img_path = joinpath(img_folder, img_file)
    try
        result = Temp_detect(img_path)
        push!(detected, img_file)
        dest_path = joinpath(detected_folder, img_file)
        mv(img_path, dest_path, force=true)
        println("✓ DETECTED: $img_file → moved to $detected_folder/")
    catch e
        push!(not_detected, img_file)
        dest_path = joinpath(not_detected_folder, img_file)
        mv(img_path, dest_path, force=true)
        println("✗ NOT DETECTED: $img_file → moved to $not_detected_folder/")
    end
end

# println()
# println("=" * 60)
# println("Summary:")
# println("  Detected:     $(length(detected))")
# println("  Not detected: $(length(not_detected))")
# println("=" * 60)

# 5 clicks wordt gedetecteerd