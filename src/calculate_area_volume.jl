

# Compute polygon area using shoelace formula
function polygon_area(coords::Vector{Tuple{Float64,Float64}})
    n = length(coords)
    area = 0.0
    for i in 1:n
        x1, y1 = coords[i]
        x2, y2 = coords[mod(i, n) + 1]  # wrap around
        area += (x1 * y2) - (x2 * y1)
    end
    return abs(area) / 2.0
end

# Compute axisymmetric volume from full outline

function axisymmetric_volume(coords::Vector{Tuple{Float64,Float64}})
    # Find axis of symmetry
    xs = [p[1] for p in coords]
    x_axis = (minimum(xs) + maximum(xs)) / 2.0

    # Group points by y and find max radius at each height
    y_to_radii = Dict{Float64, Vector{Float64}}()
    for (x, y) in coords
        r = abs(x - x_axis)
        if haskey(y_to_radii, y)
            push!(y_to_radii[y], r)
        else
            y_to_radii[y] = [r]
        end
    end

    # Create profile: (max radius, y)
    profile = [(maximum(radii), y) for (y, radii) in y_to_radii]
    sorted_profile = sort(profile, by = x -> x[2])

    # Integrate using trapezoidal rule for πr²
    volume = 0.0
    for i in 1:(length(sorted_profile)-1)
        r1, y1 = sorted_profile[i]
        r2, y2 = sorted_profile[i+1]
        avg_r2 = (r1^2 + r2^2) / 2.0
        dy = y2 - y1
        volume += π * avg_r2 * dy
    end
    return volume
end



# Function to compute axisymmetric volume for right-side profile (axis at x = 0)
function axisymmetric_volume_right(coords::Vector{Tuple{Float64,Float64}})
    # Sort points by y (height)
    sorted_coords = sort(coords, by = p -> p[2])

    # println("Sorted profile (r, y):")
    for (r, y) in sorted_coords
        # println("r = $r, y = $y")
    end

    # Integrate using trapezoidal rule for πr²
    volume = 0.0
    for i in 1:(length(sorted_coords)-1)
        r1, y1 = sorted_coords[i]
        r2, y2 = sorted_coords[i+1]
        avg_r2 = (r1^2 + r2^2) / 2.0
        dy = y2 - y1
        segment_volume = π * avg_r2 * dy
        # println("Segment $i: r1=$r1, r2=$r2, dy=$dy, avg_r²=$avg_r2, contribution=$segment_volume")
        volume += segment_volume
    end

    # println("Total Volume = $volume")
    return volume
end

function wo_fn(V, γ; g=9.81e3, Δρ=1e-3, Dn=0.9)

    wo = ( Δρ * g * V ) / (π * γ * Dn)

    return wo
end

# # Example usage:
# coords = [(1.0,0.0), (0.8,0.5), (0.5,1.0), (0.0,1.5)]  # right-side profile
# axisymmetric_volume_right(coords)


# # Example usage:
# coords = [(2.0,0.0), (2.0,1.0)]
# println("Axisymmetric Volume: ", axisymmetric_volume_right(coords))


