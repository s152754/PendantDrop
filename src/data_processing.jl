"""
Data processing functions for 
- detecting and extracting droplet edge coordinates (image coordinates)
- trimming the needle coordinates (image coordinates)
- creating a spline of the droplet edge coordinates (image coordinates)
- normalising the spline for use in the model 
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




function load_file(chosen_file)
    """
    Load an image from a chosen file or a default URL if no file is provided.

    Arguments:
    - chosen_file: The file to load, or `nothing`.

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
    """
    RectangleSelector(img, display_width=600)

    Interactive selector for defining a rectangular region of interest (ROI) on an image.

    Arguments:
    - img: The image on which to select the ROI.
    - display_width: The width at which to display the image for selection (default: 600 pixels).
    """
    
    img
    display_width::Int
end




RectangleSelector(img; display_width=600) = RectangleSelector(img, display_width)
"""
    Constructor for RectangleSelector.

    Arguments:
    - img: The image on which to select the ROI.
    - display_width: The width at which to display the image for selection (default: 600 pixels).

    Returns:
    - RectangleSelector: An instance of the RectangleSelector struct.
    """





    function img_to_dataURI(img)
    """
    Convert an image to a base64-encoded data URI for embedding in HTML.

    Arguments:
    - img: The image to convert.

    Returns:
    - String: A base64-encoded data URI representing the image.
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
    This function generates an interactive HTML canvas that allows the user to select a rectangular 
    region of interest (ROI) on the provided image. The selected ROI coordinates are stored in 
    the `value` property of the parent div element and can be accessed for further processing.

    Arguments:
    - io: The I/O stream to write to.
    - m: The MIME type for the output.
    - rs: The RectangleSelector instance to display.
    """

    # Get the dimensions of the image and calculate the display height based on the specified display width
    h, w = size(rs.img, 1), size(rs.img, 2)   # Images.jl: (rows=height, cols=width)
    scale = rs.display_width / w
    disp_w = rs.display_width
    disp_h = round(Int, h * scale)
    # Convert the image to a base64-encoded data URI for embedding in HTML
    data_uri = img_to_dataURI(rs.img)

    # Generate the HTML and JavaScript for the interactive canvas
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

function extract_roi(img, roi_raw, display_width=600)
    """
    Extract the region of interest (ROI) from the original image based on the selected rectangle coordinates.

    Arguments:
    - img: The original image from which to extract the ROI.
    - roi_raw: The raw coordinates of the selected rectangle (x0, y0, w0, h0).
    - display_width: The width at which the image was displayed for selection.

    Returns:
    - roi_img: The extracted region of interest image.
    - orig_x_min: The minimum x-coordinate of the ROI in the original image.
    - orig_y_min: The minimum y-coordinate of the ROI in the original image.
    """

	if roi_raw === nothing || length(roi_raw) != 4
		error("Draw a rectangle and click Submit to continue.")
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
			error("Selected rectangle is empty — try again.")
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
    Build a binary mask of the droplet region by combining the smoothed image and gradient magnitude.

    Arguments:
    - gradient_magnitude: The normalized gradient magnitude of the smoothed image.
    - smoothed_img: The smoothed grayscale image.
    - weight: The weight factor for combining the two signals (0 <= weight <= 1).

    Returns:
    - binary_mask: A binary mask indicating the droplet region.
    - labeled_mask: A labeled mask of connected components in the binary mask.
    - component_length: The lengths of the connected components in the labeled mask.
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





function extract_boundary_coordinates(labeled_mask; best_label=1)
    """
    Extract the boundary coordinates of the best connected component in the labeled mask.

    Arguments:
    - labeled_mask: A labeled mask of connected components.
    - best_label: The label of the connected component to extract (default: 1).

    Returns:
    - droplet_mask: A binary mask of the selected droplet region.
    - boundary_mask: A binary mask of the boundary pixels of the droplet region.     
    """

    # Select the connected component corresponding to the droplet 
    # droplet_mask is a boolean matrix 
    # true = droplet, false = background
    droplet_mask = Bool.(labeled_mask .== best_label)

    # Fill interior holes in the droplet mask to ensure a solid region
    h, w = size(droplet_mask)

    # Start from all border background pixels 
    outside_mask = falses(h, w)
    stack = CartesianIndex[]

    	for i in 1:h
	    if !droplet_mask[i, 1]
	        push!(stack, CartesianIndex(i, 1))
	    end
	    if !droplet_mask[i, w]
	        push!(stack, CartesianIndex(i, w))
	    end
	end
	
	for j in 1:w
	    if !droplet_mask[1, j]
	        push!(stack, CartesianIndex(1, j))
	    end
	    if !droplet_mask[h, j]
	        push!(stack, CartesianIndex(h, j))
	    end
	end

    # Flood fill from the border background pixels to mark all outside pixels
	while !isempty(stack)
	    c = pop!(stack)
	    i, j = Tuple(c)
	
	    if i < 1 || i > h || j < 1 || j > w
	        continue
	    end
	    if droplet_mask[i, j] || outside_mask[i, j]
	        continue
	    end
	
	    outside_mask[i, j] = true
	
	    push!(stack, CartesianIndex(i+1, j))
	    push!(stack, CartesianIndex(i-1, j))
	    push!(stack, CartesianIndex(i, j+1))
	    push!(stack, CartesianIndex(i, j-1))
	end
    
    # Remaining background pixels that are not marked as outside are considered holes inside the droplet
    holes = .!droplet_mask .& .!outside_mask

    # Combine the droplet mask with the holes to ensure a solid droplet region
    droplet_mask = droplet_mask .| holes

    # Build a contour mask from pixels adjacent to the background 
    background_mask = .!droplet_mask
    boundary_mask = droplet_mask .& (
        circshift(background_mask, ( 1,  0)) .|
        circshift(background_mask, (-1,  0)) .|
        circshift(background_mask, ( 0,  1)) .|
        circshift(background_mask, ( 0, -1))
    )

    return droplet_mask, boundary_mask

end





function trim_needle(s::Vector{<:Real}, n::Vector{<:Real}; min_points::Int=97, tol::Float64=0.5, max_consecutive_outliers::Int=55)
    """
    Trim the needle coordinates based on a linear regression model and outlier detection.

    Arguments:
    - s::Vector{<:Real}: The x-coordinates of the needle.
    - n::Vector{<:Real}: The y-coordinates of the needle.
    - min_points::Int: Minimum number of points to consider for regression (default: 90).
    - tol::Float64: Tolerance for outlier detection (default: 0.5).
    - max_consecutive_outliers::Int: Maximum number of consecutive outliers allowed (default: 55).

    Returns:
    - Int: The index of the last valid point in the trimmed needle coordinates.
    """
	#z5 = 165
	#z13: min_points::Int=70, tol::Float64=0.5, max_consecutive_outliers::Int=55)
	#z13_tilt15: min_points::Int=80, tol::Float64=0.9, max_consecutive_outliers::Int=65)

    s = Float64.(s)
    n = Float64.(n)

    if length(s) < min_points
        return 0
    end

    coeffs = hcat(s[1:min_points], ones(min_points)) \ n[1:min_points]
    cutoff = min_points
    consecutive_outliers = 0

    for i in (min_points+1):length(s)
        pred = coeffs[1] * s[i] + coeffs[2]
        if abs(n[i] - pred) > tol
            consecutive_outliers += 1
            if consecutive_outliers >= max_consecutive_outliers
                return max(cutoff - max_consecutive_outliers + 1, 0)
            end
        else
            consecutive_outliers = 0
            cutoff = i
            win = max(1, i - min_points + 1):i
            coeffs = hcat(s[win], ones(length(win))) \ n[win]
        end
    end
    
    return cutoff

end




function fit_spline(boundary_mask, orig_x_min, orig_y_min)
    """
    Fit a smoothing spline to the boundary coordinates of a droplet.

    Arguments:
    - boundary_mask: A binary mask of the droplet boundary.
    - orig_x_min: The minimum x-coordinate of the ROI in the original image.
    - orig_y_min: The minimum y-coordinate of the ROI in the original image.

    Returns:
    - xs_open: The x-coordinates of the open boundary points.
    - ys_open: The y-coordinates of the open boundary points.
    - x_contour: The x-coordinates of the contour points.
    - y_contour: The y-coordinates of the contour points.
    - x_spline_pts: The x-coordinates of the fitted spline points.
    - y_spline_pts: The y-coordinates of the fitted spline points.
    - x_spline_norm_pts: The x-coordinates of the normalised spline points.
    - y_spline_norm_pts: The y-coordinates of the normalised spline points.
    - x_spline_norm_pts_lhs: The x-coordinates of the left half of the normalised spline points.
    - y_spline_norm_pts_lhs: The y-coordinates of the left half of the normalised spline points.
    - x_spline_norm_pts_rhs: The x-coordinates of the right half of the normalised spline points.
    - y_spline_norm_pts_rhs: The y-coordinates of the right half of the normalised spline points.   
    """

	# TODO: kwargs to put in function later: row_tol, dedup_tol, λ, n_samples
	
    # STEP 1: get the FULL boundary in order

		# transpose mask so axes are [x, y], which matches what Contour.jl expects
    mask_xy = permutedims(Float64.(boundary_mask))
    cl = Contour.contours(1:size(mask_xy, 1), 1:size(mask_xy, 2), mask_xy, [0.5])
    line = first(Contour.lines(first(Contour.levels(cl))))
    xs_c, ys_c = Contour.coordinates(line)
    xs_c = xs_c .+ orig_x_min .- 1
    ys_c = ys_c .+ orig_y_min .- 1

		# marching squares closes the loop (first == last) 
    if (xs_c[1], ys_c[1]) == (xs_c[end], ys_c[end])
        pop!(xs_c); pop!(ys_c)
    end

	
    # STEP 2: make sure there's no horizontal line at the top
	row_tol = 1.5 # tolerance for the crop row(s) -- TODO: this should be a function input
    y_top = minimum(ys_c)
    on_top = ys_c .<= y_top + row_tol    
    keep_pts = .!on_top
	
	    # rotate the circular array so it starts at a kept point, then take the first contiguous run of kept points -> opens the loop
    start_idx = findfirst(keep_pts)
    keep_r = circshift(keep_pts, 1 - start_idx)
    xs_r    = circshift(xs_c, 1 - start_idx)
    ys_r    = circshift(ys_c, 1 - start_idx)
    run_end = something(findfirst(!, keep_r), length(keep_r) + 1) - 1

    xs_open = xs_r[1:run_end]
    ys_open = ys_r[1:run_end]
    pts_open = collect(zip(xs_open, ys_open))


    # STEP 3: trim needle stems from EACH END of the ordered sequence independently
    trim_cx, trim_cy = mean(xs_open), mean(ys_open)
    trim_dx, trim_dy = xs_open .- trim_cx, ys_open .- trim_cy
    trim_cov = [sum(trim_dx.^2) sum(trim_dx.*trim_dy); sum(trim_dx.*trim_dy) sum(trim_dy.^2)] ./ length(xs_open)
    trim_vals, trim_vecs = eigen(trim_cov)
    trim_major = trim_vecs[:, argmax(trim_vals)] ./ norm(trim_vecs[:, argmax(trim_vals)])
    trim_perp  = [-trim_major[2], trim_major[1]]
    trim_n = trim_dx .* trim_perp[1] .+ trim_dy .* trim_perp[2]   # "straightness" coordinate

    d_open = sqrt.(diff(xs_open).^2 .+ diff(ys_open).^2)
    s_fwd = vcat(0.0, cumsum(d_open))          # arc length from the start
    s_bwd = vcat(0.0, cumsum(reverse(d_open))) # arc length from the end

    cut_start = trim_needle(s_fwd, trim_n)
    cut_end   = trim_needle(s_bwd, reverse(trim_n))

    lo = cut_start + 1
    hi = length(pts_open) - cut_end
    contour_pts = pts_open[lo:hi]

    x_contour = first.(contour_pts)
    y_contour = last.(contour_pts)

	
	# STEP 4: fit spline to contour
	
		# Dedupe: remove near-zero-distance consecutive points
		# (Dierckx needs strictly increasing arc-length knots; jagged pixel boundaries can create these)
    dedup_tol = 1e-6 # TODO: make this a function input
    step_len = hypot.(diff(x_contour), diff(y_contour))
    keep_spline = vcat(true, step_len .> dedup_tol)

    xs_dedup = x_contour[keep_spline]
    ys_dedup = y_contour[keep_spline]

    	# build cumulative arc-length parameterization
    ds = sqrt.(diff(xs_dedup).^2 .+ diff(ys_dedup).^2)
    t  = vcat(0.0, cumsum(ds))

    	# fit smoothing splines (s = λ), not exact interpolation
    λ = 5   # tune: larger = smoother/less jagged, but can round off the tip TODO: make this a function input
    x_spline = Spline1D(t, xs_dedup; s=λ)
    y_spline = Spline1D(t, ys_dedup; s=λ)

    	# sample along open arc only
    n_samples = 1200 # TODO: make this a function input
    s_sample = range(0.0, t[end], length=n_samples)
    x_spline_pts = x_spline.(s_sample)
    y_spline_pts = y_spline.(s_sample)

	
    # STEP 5: Split spline into left / right halves
	    # 1. rotate so p1 and p2 share the same y
	    # 2. translate so their midpoint sits at x=0, y=0
	    # 3. split left/right at x=0
    p1 = (x_spline_pts[1],   y_spline_pts[1])
    p2 = (x_spline_pts[end], y_spline_pts[end])

    	# angle needed to make the p1->p2 vector horizontal
    ang = atan(p2[2]-p1[2], p2[1]-p1[1])
    cA, sA = cos(-ang), sin(-ang)

    rotate_pt(x, y) = (x*cA - y*sA, x*sA + y*cA)

    rot_pts = rotate_pt.(x_spline_pts, y_spline_pts)
    xr = first.(rot_pts)
    yr = last.(rot_pts)

    mid_x = (xr[1] + xr[end]) / 2
    mid_y = (yr[1] + yr[end]) / 2

    x_spline_norm_pts = xr .- mid_x
    y_spline_norm_pts = yr .- mid_y

		# check droplet orientation is correct
    if mean(y_spline_norm_pts) < 0
        x_spline_norm_pts = .-x_spline_norm_pts
        y_spline_norm_pts = .-y_spline_norm_pts
    end
    	# split at x = 0
    mask_lhs = x_spline_norm_pts .<= 0
    mask_rhs = .!mask_lhs

    x_spline_norm_pts_lhs, y_spline_norm_pts_lhs = x_spline_norm_pts[mask_lhs], y_spline_norm_pts[mask_lhs]
    x_spline_norm_pts_rhs, y_spline_norm_pts_rhs = x_spline_norm_pts[mask_rhs], y_spline_norm_pts[mask_rhs]

    # TODO: make a reverse coordinate translation for the RHS

    return xs_open, ys_open, x_contour, y_contour, x_spline_pts, y_spline_pts, x_spline_norm_pts, y_spline_norm_pts, x_spline_norm_pts_lhs, y_spline_norm_pts_lhs, x_spline_norm_pts_rhs, y_spline_norm_pts_rhs
end