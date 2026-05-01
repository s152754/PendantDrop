"""
Draw a radius line from the circle center (cx, cy) to the perimeter at angle θ.

Arguments:
- plt:   a Plots.jl plot object to mutate
- cx, cy: center coordinates
Keyword arguments:
- r::Real=1.0        radius (length of the line)
- θ::Real=0.0        angle in radians (0 .. 2π) measured from +x axis, CCW
- color=:black
- lw::Real=1.5       line width
- linestyle=:solid   :solid, :dash, :dot, ...
- cap=:butt          line cap style (:butt, :round, :square) if supported
"""
function radius_line!(plt, cx, cy, r, θ;
                      color=:black, lw::Real=1.5,
                      linestyle::Symbol=:solid, cap=:butt)

    # Normalize θ to [0, 2π) to be safe
    θn = mod(θ, 2π)

    # Endpoint on the circle
    x_end = cx + r * cos(θn)
    y_end = cy + r * sin(θn)

    # Draw the line from center to endpoint
    plot!(plt, [cx, x_end], [cy, y_end];
          color=color, lw=lw, linestyle=linestyle, label="")
    return plt
end


"""
Draw a dashed circle on an existing plot `plt`.

Arguments:
- plt:        a Plots.jl plot object to mutate
- cx, cy:     center coordinates of the circle
Keyword arguments:
- r::Real=1.0       radius
- start_angle=0.0   start angle (radians)
- end_angle=2π      end angle (radians)
- n::Int=200        number of points
- color=:black
- lw::Real=1      line width
"""
function dashed_circle!(plt, cx, cy, r;
                        start_angle::Real=0.0,
                        end_angle::Real=2π,
                        n::Int=200,
                        color=:black,
                        lw::Real=1)
    θ = range(start_angle, end_angle; length=n)
    xs = cx .+ r .* cos.(θ)
    ys = cy .+ r .* sin.(θ)
    plot!(plt, xs, ys; color=color, linestyle=:dash, lw=lw, label="")
    return plt
end


function dashed_arc!(plt, cx, cy; r=0.15, θ1=0.0, θ2=pi/2, n=60, color=:black)
    θ = range(θ1, θ2; length=n)
    xs = cx .+ r .* cos.(θ)
    ys = cy .+ r .* sin.(θ)
    plot!(plt, xs, ys, color=color, linestyle=:dash, lw=1, label="")
end

"""
Draw a dashed ellipse (tilted circle) on a given plot.

Arguments:
- plt: the plot object
- cx, cy: center coordinates
- r: radius of the original circle
- scale_y: scaling factor for y-axis (e.g., 0.6 for flattening)
- angle: rotation in radians (e.g., π/6 for tilt)
- n: number of points for smoothness
"""
function dashed_ellipse!(plt, cx, cy; r=1.0, scale_y=0.6, angle=0.0, start_angle=0.0, end_angle=2π, n=200, color=:black)
    θ = range(start_angle, end_angle; length=n)
    xs = r .* cos.(θ)
    ys = r .* sin.(θ) .* scale_y  # flatten vertically

    # Apply rotation
    rot_x = xs .* cos(angle) .- ys .* sin(angle)
    rot_y = xs .* sin(angle) .+ ys .* cos(angle)

    # Shift to center
    xs_final = cx .+ rot_x
    ys_final = cy .+ rot_y

    plot!(plt, xs_final, ys_final, color=color, linestyle=:dash, lw=1, label="")
end



"""
Draw gray x-axis and y-axis lines, then add short arrows with labels 'x' and 'y'.

Arguments:
- plt: plot object
- xlims, ylims: axis limits
- arrow_len: length of the short arrows
- color_axes: color for axes lines
- color_arrows: color for arrows
"""
function add_axes_with_arrows!(plt; xlims=(-1, 1), ylims=(-1, 1),
                               arrow_len=0.2, color_axes=:black, color_arrows=:black)

    # Gray axes lines
    plot!(plt, [xlims[1], xlims[2]], [0,0], color=color_axes, lw=1, label="", linestyle=:dash)  # x-axis
    plot!(plt, [0, 0], [ylims[1], ylims[2]*1.2], color=color_axes, lw=1, label="", linestyle=:dash)  # y-axis
end

function add_axes_with_arrows_short!(plt; xlims=(-1, 1), ylims=(-1, 1), arrow_len=0.2,
                               color=:black)

    # Short arrow on x-axis
    plot!(plt, [0, arrow_len], [ylims[2], ylims[2]], color=:black, lw=1, label="", arrow=:closed)
    annotate!(plt, arrow_len, ylims[2]-0.18, text(L"r", :black))

    # Short arrow on y-axis
    plot!(plt, [0, 0], [ylims[2], ylims[2]+arrow_len], color=:black, lw=1, label="", arrow=:closed)
    annotate!(plt, -0.1, ylims[2]+arrow_len + 0.05, text(L"z", :black))
end



function tangent_to_y0(curve::Matrix{Float64}, idx::Int)
    # Extract point
    x0, y0 = curve[idx, :]

    # Compute slope using central difference
    m = (curve[idx+1, 2] - curve[idx-1, 2]) / (curve[idx+1, 1] - curve[idx-1, 1])

    # Compute intersection with y = 0
    x_end = x0 - y0 / m

    # Return start and end points
    return [x0, x_end], [y0, 0.0]
end



"""
    initial_y(theta, endpoint; angle_from=:x, degrees=false)

Compute the y-coordinate `y0` of the initial point on the y-axis (x=0) such that
a straight line at angle `theta` reaches the given `endpoint` (x_end, y_end).

Arguments:
- `theta`::Real             Angle of the line.
- `endpoint`::Tuple{<:Real,<:Real}  (x_end, y_end) point on the curve the line must pass through.

Keyword arguments:
- `angle_from`::Symbol      :x  -> theta measured from x-axis (default, slope = tan(theta))
                            :y  -> theta measured from y-axis (slope = cot(theta))
- `degrees`::Bool           If true, interpret `theta` in degrees (converted to radians).

Returns:
- `y0`::Float64             y-coordinate at x=0 where the line starts on the y-axis.

Notes / Edge cases:
- If `angle_from=:x` and `theta ≈ π/2` (or 90°), the line is vertical and does NOT intersect the y-axis
  (undefined `y0`); an error is thrown.
- If `angle_from=:y` and `theta ≈ 0` (or 0°), the line is vertical and does NOT intersect the y-axis
  (undefined `y0`); an error is thrown.
"""
function initial_y(theta::Real, endpoint::Tuple{<:Real,<:Real}; angle_from::Symbol=:x, degrees::Bool=false)
    x_end, y_end = float(endpoint[1]), float(endpoint[2])

    # Convert to radians if needed
    θ = degrees ? (theta * π / 180) : float(theta)

    # Compute slope m based on angle reference
    m = if angle_from == :x
        # slope = tan(theta)
        # Guard against vertical line (tan -> Inf)
        abs(cos(θ)) < 1e-12 && error("Line is vertical for θ ≈ π/2; it does not intersect the y-axis.")
        tan(θ)
    elseif angle_from == :y
        # slope = cot(theta) = tan(π/2 - θ)
        abs(sin(θ)) < 1e-12 && error("Line is vertical for θ ≈ 0; it does not intersect the y-axis.")
        tan(π/2 - θ)
    else
        error("angle_from must be :x or :y")
    end

    # y0 from y_end = y0 + m * x_end
    y0 = y_end - m * x_end
    return y0
end


"""
    reorder_droplet(x::Vector{Int}, y::Vector{Int}) -> x_new::Vector{Int}, y_new::Vector{Int}

Given pixel coordinates `x` and `y` (same length), this function:
1) Sorts points by x (left to right) and splits into left and right halves.
2) Sorts the left edge by y descending (top → bottom).
3) Sorts the right edge by y ascending (bottom → top).
4) Concatenates the two edges and returns reordered x and y as integers.

This is useful to order boundary points of a roughly “droplet-shaped” contour.
"""
function reorder_dropletpx(x::Vector{Int}, y::Vector{Int})
    @assert length(x) == length(y) "x and y must have the same length"

    # Combine into tuples of (x, y)
    points = collect(zip(x, y))

    # Sort by x to split left/right
    sorted_points = sort(points; by = p -> p[1])

    # Split roughly in half (left half gets the smaller x's)
    n = length(sorted_points)
    mid = n ÷ 2                       # floor division
    left_edge  = view(sorted_points, 1:mid)
    right_edge = view(sorted_points, mid+1:n)

    # Left edge: sort by y descending (top → bottom if y increases downward)
    left_sorted  = sort(left_edge;  by = p -> p[2], rev = true)

    # Right edge: sort by y ascending (bottom → top)
    right_sorted = sort(right_edge; by = p -> p[2], rev = false)

    # Concatenate
    reordered = vcat(right_sorted, left_sorted)

    # Unzip back to integer vectors
    x_new = [p[1] for p in reordered]
    y_new = [p[2] for p in reordered]

    return x_new, y_new
end


function reorder_droplet(x::Vector{Float64}, y::Vector{Float64})
    # Combine into tuples
    points = collect(zip(x, y))

    # Sort by x to split left/right
    sorted_points = sort(points, by = p -> p[1])

    # Split roughly in half
    mid = Int(length(points) ÷ 2)
    left_edge = sorted_points[1:mid]
    right_edge = sorted_points[mid+1:end]

    # Sort left edge by y descending (top to bottom)
    left_edge = sort(left_edge, by = p -> p[2], rev = true)

    # Sort right edge by y ascending (bottom to top)
    right_edge = sort(right_edge, by = p -> p[2])

    # Concatenate
    reordered = vcat(left_edge, right_edge)

    # Return separate x and y
    x_new = [p[1] for p in reordered]
    y_new = [p[2] for p in reordered]

    return x_new, y_new
end


function mycircle!(plt, cx, cy, r;
                        start_angle::Real=0.0,
                        end_angle::Real=2π,
                        n::Int=200,
                        color=:black,
                        lw::Real=1.5,
                        )
    θ = range(start_angle, end_angle; length=n)
    xs = cx .+ r .* cos.(θ)
    ys = cy .+ r .* sin.(θ)
    plot!(plt, xs, ys; color=color, lw=lw, label=latexstring("R_\\mathrm{apex} = $(r) \\ \\mathrm{mm}")  )
    return plt
end

function QoIppd_plotting(Nsamples, samplesppd, QoIppd, flatchain)
    for i in 1:Nsamples
        samplesppd[i] = QoIppd[i] .- QoIppd[i][end,2]
    end


    meanθs = mean(flatchain, dims=2)
    rmmean, zmmean = modelfn(vec(meanθs))
    meanshape = hcat([rmmean, zmmean]...)

    alignedpoint = align_to_reference(meanshape, samplesppd)
    normalsref = compute_normals_2d(meanshape)
    lower_curve, upper_curve = credible_band(meanshape, normalsref, alignedpoint)

    x_mean = meanshape[:, 1]
    y_mean = meanshape[:, 2]

    x_upper = upper_curve[:, 1]
    y_upper = upper_curve[:, 2]

    x_lower = lower_curve[:, 1]
    y_lower = lower_curve[:, 2]

    return x_mean, y_mean, x_upper, y_upper, x_lower, y_lower

end
