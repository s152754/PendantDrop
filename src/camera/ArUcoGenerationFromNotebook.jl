using PyCall
using FileIO
using Images

function generate_aruco_marker(marker_id::Int, marker_size_cm::Float64, dpi::Int; 
                               dictionary::String="DICT_4X4_50", output_dir::String=".")
    """
    Generate an ArUco marker and save as PNG.
    
    Args:
        marker_id: ID of the marker (must be valid for chosen dictionary)
        marker_size_cm: Physical size of marker in centimeters (side length)
        dpi: Print resolution in dots per inch
        dictionary: ArUco dictionary to use (default: DICT_4X4_50)
        output_dir: Directory to save the marker image
    
    Returns:
        Path to saved image file
    """
    
    cv2 = pyimport("cv2")
    aruco = cv2[:aruco]
    
    # Convert dictionary string to OpenCV constant
    dict_name = Symbol(dictionary)
    dict_constant = getfield(aruco, dict_name)
    
    # Get the dictionary
    aruco_dict = aruco.getPredefinedDictionary(dict_constant)
    
    # Compute pixel size for desired real-world size
    # 1 inch = 2.54 cm, so pixels = (cm / 2.54) * dpi
    marker_size_px = Int(round((marker_size_cm / 2.54) * dpi))
    
    # Generate the ArUco marker
    marker_img = aruco.generateImageMarker(aruco_dict, marker_id, marker_size_px)
    
    # Create output directory if it doesn't exist
    if !isdir(output_dir)
        mkpath(output_dir)
    end
    
    # Save the marker image
    filename = joinpath(output_dir, "aruco_marker_$(marker_id)_$(marker_size_cm)cm.png")
    save(filename, marker_img)
    
    # Display information
    println("✓ Marker generated successfully")
    println("  ID: $marker_id")
    println("  Size: $marker_size_cm cm × $marker_size_cm cm")
    println("  Pixels: $marker_size_px px at $dpi dpi")
    println("  Dictionary: $dictionary")
    println("  Saved to: $filename")
    
    return filename
end

function print_available_dictionaries()
    """Print available ArUco dictionaries"""
    cv2 = pyimport("cv2")
    aruco = cv2[:aruco]
    
    # Get all dictionary constants
    dicts = [
        "DICT_4X4_50",
        "DICT_4X4_100",
        "DICT_4X4_250",
        "DICT_4X4_1000",
        "DICT_5X5_50",
        "DICT_5X5_100",
        "DICT_5X5_250",
        "DICT_5X5_1000",
        "DICT_6X6_50",
        "DICT_6X6_100",
        "DICT_6X6_250",
        "DICT_6X6_1000",
        "DICT_7X7_50",
        "DICT_7X7_100",
        "DICT_7X7_250",
        "DICT_7X7_1000",
        "DICT_ARUCO_ORIGINAL",
    ]
    
    println("Available ArUco Dictionaries:")
    for dict_name in dicts
        try
            dict_constant = getfield(aruco, Symbol(dict_name))
            aruco_dict = aruco.getPredefinedDictionary(dict_constant)
            n_markers = aruco_dict.bytesList .* 8  # Number of bits per marker
            println("  - $dict_name")
        catch
        end
    end
end

function print_usage()
    println("Usage: julia ArUcoGenerationFromNotebook.jl [marker_id] [size_cm] [dpi] [output_dir]")
    println()
    println("Arguments:")
    println("  marker_id    : ID of marker (default: 7)")
    println("  size_cm      : Physical size in centimeters (default: 0.5)")
    println("  dpi          : Print resolution in DPI (default: 300)")
    println("  output_dir   : Directory to save marker (default: current directory)")
    println()
    println("Examples:")
    println("  julia ArUcoGenerationFromNotebook.jl")
    println("  julia ArUcoGenerationFromNotebook.jl 5 2.0 300")
    println("  julia ArUcoGenerationFromNotebook.jl 10 5.0 300 ./markers")
    println()
    println("Note: You can also use https://chev.me/arucogen/ for marker generation")
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Default values
    default_marker_id = 7
    default_marker_size_cm = 0.5
    default_dpi = 300
    default_output_dir = "."
    default_dictionary = "DICT_4X4_50"
    
    if length(ARGS) > 0 && ARGS[1] in ["-h", "--help"]
        print_usage()
        exit()
    end
    
    if length(ARGS) > 0 && ARGS[1] in ["--list-dicts", "--dictionaries"]
        print_available_dictionaries()
        exit()
    end
    
    # Parse command-line arguments
    marker_id = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : default_marker_id
    marker_size_cm = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : default_marker_size_cm
    dpi = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : default_dpi
    output_dir = length(ARGS) >= 4 ? ARGS[4] : default_output_dir
    
    println("ArUco Marker Generation")
    println("=" * 50)
    println("Parameters:")
    println("  Marker ID: $marker_id")
    println("  Size: $marker_size_cm cm")
    println("  DPI: $dpi")
    println("  Dictionary: $default_dictionary")
    println("  Output dir: $output_dir")
    println()
    
    try
        generate_aruco_marker(marker_id, marker_size_cm, dpi; 
                            dictionary=default_dictionary, 
                            output_dir=output_dir)
    catch e
        println("✗ Error generating marker:")
        println("  $(e)")
        exit(1)
    end
end
