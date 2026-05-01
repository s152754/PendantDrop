using JLD2, Plots

case = "14nov_thres5"

# compare orthographic to perspective
Lp = load("data/dropletpx_$(case).jld2")
rp = Lp["rpx"]
zp = -Lp["Zpx"]

function sortpx_fn(rp, zp)
    # pair r and z into point tuples
    points = collect(zip(rp, zp))

    # choose a splitting threshold (fall back to last element if too short)
    thresh = length(rp) >= 2 ? rp[end-1] : rp[end]

    # filter properly: predicate first, collection second
    left_side  = filter(p -> p[1] <= thresh, points)
    right_side = filter(p -> p[1]  > thresh, points)

    # sort each side by z (y) descending
    left_sorted  = sort(left_side,  by = p -> p[2], rev = true)
    right_sorted = sort(right_side, by = p -> p[2], rev = true)

    # build closed polygon: left (downwards), right (upwards), close loop
    polygon = vcat(left_sorted, reverse(right_sorted))
    push!(polygon, polygon[1])

    poly_x = first.(polygon)
    poly_y = last.(polygon)

    return poly_x, poly_y
end


rps, zps = sortpx_fn(rp, zp)

figcomp = plot(title="Orthographic vs Perspective", xlabel="r (mm)", ylabel="z (mm)", legend=:topright, aspect_ratio=:equal, grid=true)
scatter!(figcomp, rps, zps; markershape=:diamond, label="perspective")
display(figcomp)

# create Surface
rp_closed = [rps; rps[1]]
zp_closed = [zps; zps[1]]

figclose = plot(rp_closed, zp_closed, fill=(0, :gray); title="Closed Surface from Perspective Data", xlabel="r (mm)", ylabel="z (mm)", legend=false, aspect_ratio=:equal, grid=true)
plot!(figclose, axis = false,    # removes axes & ticks
     grid = false,    # removes grid
     framestyle = :none,  # removes the frame
     legend = false,  # remove legend
     title = "" )      # remove title
display(figclose)

savefig(figclose, "data/filled/$(case).png")