"""
Data processing functions for detecting and extracting droplet edge coordinates from uploading image
"""

function test(x::Int)
    """
    Test function for the PendantDropModule.
    
    Arguments:
    - x::Int: An integer input.

    Returns:
    - Int: The input integer incremented by 3.
    """

    return x + 3

end

function load_file(chosen_file::Union{Nothing, Dict{String, Any}})
    """
    Load an image from a chosen file or a default URL if no file is provided.

    Arguments:
    - chosen_file::Union{Nothing, Dict{String, Any}}: A dictionary containing the chosen file's data or `nothing`.

    Returns:
    - img: The loaded image.
    """
	if chosen_file !== nothing
		img_bytes = chosen_file["data"] 
		img = load(IOBuffer(img_bytes))
	else 
		img_url = "https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcSpSdtJ3NdDYxIp5-Abe9PifEG_y-fDQP-QchWxLoLXMMpSfVtgRy0FKYS1GmWTjxx_QRs&usqp=CAU"
		img_file = download(img_url)
		img = load(img_file)
	end

    return img

end

struct RectangleSelector
    img
    display_width::Int
end

RectangleSelector(img; display_width=600) = RectangleSelector(img, display_width)

function img_to_datauri(img)
    """

    """
    io = IOBuffer()
    # writing via a temp file is the most reliable path across FileIO/ImageIO versions
    tmp = tempname() * ".png"
    save(tmp, img)
    b64 = base64encode(read(tmp))
    rm(tmp, force=true)
    return "data:image/png;base64,$(b64)"
end

function Base.show(io::IO, m::MIME"text/html", rs::RectangleSelector)
    """

    """

    h, w = size(rs.img, 1), size(rs.img, 2)   # Images.jl: (rows=height, cols=width)
    scale = rs.display_width / w
    disp_w = rs.display_width
    disp_h = round(Int, h * scale)
    data_uri = img_to_datauri(rs.img)

    Base.show(io, m, @htl("""
    <div>
    <canvas width="$(disp_w)" height="$(disp_h)" style="cursor: crosshair; touch-action: none;"></canvas>
    <script>
    const div = currentScript.parentElement
    div.value = [0, 0, 0, 0]
    
    const canvas = div.querySelector("canvas")
    const ctx = canvas.getContext("2d")
    
    const image = new Image()
    image.src = $(data_uri)
    
    var startX = 0, startY = 0, curW = 0, curH = 0
    
    function draw() {
        ctx.clearRect(0, 0, canvas.width, canvas.height)
        ctx.drawImage(image, 0, 0, canvas.width, canvas.height)
        ctx.fillStyle = 'rgba(63, 61, 109, 0.3)'
        ctx.strokeStyle = '#ff3b3b'
        ctx.lineWidth = 2
        ctx.fillRect(startX, startY, curW, curH)
        ctx.strokeRect(startX, startY, curW, curH)
    }
    
    image.onload = () => draw()
    
    function onmove(e) {
        curW = e.offsetX - startX
        curH = e.offsetY - startY
        div.value = [startX, startY, curW, curH]
        div.dispatchEvent(new CustomEvent("input"))
        draw()
    }
    
    canvas.onpointerdown = e => {
        startX = e.offsetX
        startY = e.offsetY
        curW = 0
        curH = 0
        canvas.onpointermove = onmove
    }
    
    window.addEventListener("pointerup", e => {
        canvas.onpointermove = null
    })
    </script>
    </div>
    """))
end

function extract_roi(img, roi_raw, display_width)
    """

    """
	if roi_raw === nothing || length(roi_raw) != 4
		roi_img = img
		orig_x_min, orig_y_min = 1, 1
	else
		x0, y0, w0, h0 = roi_raw

		x_min, x_max = minmax(x0, x0 + w0)
		y_min, y_max = minmax(y0, y0 + h0)

		scale = size(img, 2) / display_width

		orig_x_min = clamp(round(Int, x_min * scale) + 1, 1, size(img, 2))
		orig_x_max = clamp(round(Int, x_max * scale),     1, size(img, 2))
		orig_y_min = clamp(round(Int, y_min * scale) + 1, 1, size(img, 1))
		orig_y_max = clamp(round(Int, y_max * scale),     1, size(img, 1))

		if orig_x_max > orig_x_min && orig_y_max > orig_y_min
			roi_img = img[orig_y_min:orig_y_max, orig_x_min:orig_x_max]
		else
			roi_img = img
		end
	end

    return roi_img, orig_x_min, orig_y_min

end

function find_max_gradients(roi_img)
    """
    Compute the gradient magnitude and smoothed image from the region of interest (ROI) image.

    Arguments:
    - roi_img: The region of interest image.

    Returns:
    - gradient_magnitude: The normalized gradient magnitude of the smoothed image.
    - smoothed_img: The smoothed grayscale image.
    """

	# Convert the ROI image to grayscale
	grayscale_img = Gray.(roi_img)

	# Smooth the image to reduce noise before edge detection
	smoothed_img = imfilter(grayscale_img, Kernel.gaussian(2.0))
	smoothed_img = clamp01.(smoothed_img) # normalise
	
	# Compute intensity gradients in the x and y directions 
	gradient_x, gradient_y = imgradients(smoothed_img, KernelFactors.prewitt, "replicate")

    # Combine gradients into a single edge detection strength map
	gradient_magnitude = sqrt.(gradient_x.^2 .+ gradient_y.^2)

    # Normalise the gradient magnitude to range [0, 1]
	gradient_magnitude ./= maximum(gradient_magnitude) 

    return gradient_magnitude, smoothed_img

end

function build_binary_mask(gradient_magnitude, smoothed_img, weight)
    """

    """

    # Combine the two signals into a single score image
    score_img = weight .* (1 .- smoothed_img) .+ (1-weight) .* gradient_magnitude

    # Estimate a threshold using Otsu's method
    threshold = find_threshold(score_img, Otsu())

    # Create a binary mask from the score image and threshold
    binary_mask = score_img .> threshold

    # Idenitfy connected components in the binary mask
    labeled_mask = label_components(binary_mask)
    component_length = component_lengths(labeled_mask)

    return binary_mask, labeled_mask, component_length

end

function extract_boundary_coordinates(labeled_mask, best_label, offset; contour_thickness=1)
    """
    Extract the boundary coordinates of the best connected component in the labeled mask.

 
    """

    # Select the connected component corresponding to the droplet 
    droplet_mask = Bool.(labeled_mask .== best_label)

    # Build a contour mask from pixels adjacent to the background 
    background_mask = .!droplet_mask
    contour_mask = droplet_mask .& (
        circshift(background_mask, (1, 0)) .|
        circshift(background_mask, (-1, 0)) .|
        circshift(background_mask, (0, 1)) .|
        circshift(background_mask, (0, -1))
    )

    # Compute a signed distance transform around the contour 
    # Compute a signed distance transform around the contour
    distance_inside = distance_transform(feature_transform(.!contour_mask))
    distance_outside = distance_transform(feature_transform(contour_mask))
    signed_distance = distance_outside .- distance_inside

    # Create a shifted contour at the requested offset
    shifted_contour = abs.(signed_distance .- offset) .< contour_thickness

    # Convert to a boolean mask for downstream morphology
    binary_mask = Bool.(shifted_contour)

    # Close small gaps and isolate boundary pixels
    closed_mask = closing(binary_mask; dims=(2,))
    eroded_mask = erode(closed_mask)
    boundary_mask = closed_mask .& .!eroded_mask

    return droplet_mask, boundary_mask

end
