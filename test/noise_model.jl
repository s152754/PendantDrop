using Random
using Plots
using NearestNeighbors

# Random.seed!(42)

# --- Step 1: Generate Perfect Circle Points ---
function generate_circle(center::Tuple{Float64,Float64}, radius::Float64, n_points::Int)
    x_c, y_c = center
    xs = [x_c + radius * cos(θ) for θ in range(0, 2π, length=n_points)]
    ys = [y_c + radius * sin(θ) for θ in range(0, 2π, length=n_points)]
    return xs, ys
end

# --- Step 2: Generate Noisy Points ---
function generate_noisy_points(center::Tuple{Float64,Float64}, radius::Float64, n_points::Int; noise_std=0.3)
    x_c, y_c = center
    xs = Float64[]
    ys = Float64[]
    for i in 1:n_points
        θ = 2π * rand()
        x = x_c + radius * cos(θ)
        y = y_c + radius * sin(θ)
        δ = randn() * noise_std
        push!(xs, x + δ * cos(θ))
        push!(ys, y + δ * sin(θ))
    end
    return xs, ys
end

# --- Step 3: Build KD-Tree and Find Nearest Neighbors ---
function find_nearest_neighbors(xs_circle, ys_circle, xs_noisy, ys_noisy)
    perfect_points = hcat(xs_circle, ys_circle)'  # 2×N
    noisy_points = hcat(xs_noisy, ys_noisy)'      # 2×M
    tree = KDTree(perfect_points)
    idxs, dists = knn(tree, noisy_points, 1)      # 1 nearest neighbor per noisy point
    return getindex.(idxs, 1), dists  # Convert vector of vectors to flat vector
end

# --- Step 4: Plot Everything ---
function plot_projection(xs_circle, ys_circle, xs_noisy, ys_noisy, idxs)
    scatter(xs_circle, ys_circle, label="Perfect Circle", color=:blue, aspect_ratio=:equal)
    scatter!(xs_noisy, ys_noisy, label="Noisy Points", color=:red)
    # Highlight closest neighbors
    xs_closest = xs_circle[idxs]
    ys_closest = ys_circle[idxs]
    scatter!(xs_closest, ys_closest, label="Closest Neighbors", color=:green, markersize=6, legend=:bottomright)
    # Draw lines from noisy points to their closest neighbor
    for i in 1:length(xs_noisy)
        plot!([xs_noisy[i], xs_closest[i]], [ys_noisy[i], ys_closest[i]], color=:gray, alpha=0.5, label="")
    end
    title!("Noisy Points and Their Closest Perfect Points")
end

# --- Main Execution ---
center = (0.0, 0.0)
radius = 5.0



figsinglep = plot(xlabel="npoints", ylabel="distance", title="Distance to Closest Point vs Number of Points", xaxis=:log)

npoints = [25, 50, 100, 200, 400, 1600, 3200]
npoints = [1,10,20, 40, 80, 160, 320, 640, 1280]
seeds = [1, 4, 42, 10, 500]
colorseeds = [:red, :green, :blue, :orange, :purple]
figproji= plot(title="Noisy Points and Their Closest Perfect Points (n=25)", xlabel="X", ylabel="Y", aspect_ratio=:equal, legend=:topleft)
for j in eachindex(seeds)
    Random.seed!(seeds[j])
   

    # overkill points
    nend = npoints[end]
    xs_circleend, ys_circleend = generate_circle(center, radius, nend)
    xs_noisy1, ys_noisy1 = generate_noisy_points(center, radius, nend; noise_std=0.5)
    idxsend, distsend = find_nearest_neighbors(xs_circleend, ys_circleend, xs_noisy1, ys_noisy1)
    

    for i in eachindex(npoints)
        n = npoints[i]
        xs_noisy, ys_noisy = generate_noisy_points(center, radius, n; noise_std=0.5)
        xs_circle, ys_circle = generate_circle(center, radius, npoints[end])
        
        idxs, dists = find_nearest_neighbors(xs_circle, ys_circle, xs_noisy, ys_noisy)
        # use scalar ratio (mean of distances) instead of vector-matrix division
        scatter!(figsinglep, [n], [mean(dists)/mean(distsend)], label="", color=colorseeds[j])
        println("Seed $(seeds[j]), n=$(n): Mean Distance = $(mean(dists))")
        # plot_projection(xs_circle, ys_circle, xs_noisy, ys_noisy, idxs)
        # avg_dist = mean(dists)
        # scatter!(figsinglep, [n], [avg_dist], label="", color=:blue)


        # show circle with points
        figproj = plot()
        scatter!(figproj, xs_circle, ys_circle, label="Perfect Circle", color=:blue, aspect_ratio=:equal)
        scatter!(figproj, xs_noisy, ys_noisy, label="Noisy Points", color=colorseeds[j])
        # Highlight closest neighbors
        xs_closest = xs_circle[idxs]
        ys_closest = ys_circle[idxs]
        scatter!(figproj, xs_closest, ys_closest, label="Closest Neighbors", color=:green, markersize=6, legend=:bottomright)
        # Draw lines from noisy points to their closest neighbor
        for i in 1:length(xs_noisy)
            plot!(figproj, [xs_noisy[i], xs_closest[i]], [ys_noisy[i], ys_closest[i]], color=:gray, alpha=0.5, label="")
        end
        title!(figproj, "Noisy Points and Their Closest Perfect Points, n=$(n))")
        display(figproj)


        if i == 1
            if j == 1
            scatter!(figproji, xs_circle, ys_circle, label="Perfect Circle", color=:black, aspect_ratio=:equal)
            end
            scatter!(figproji, xs_noisy, ys_noisy, label="noisy point $(j)", color=colorseeds[j])
        end

    end
    
end


display(figproji)

display(figsinglep)
# xs_circle, ys_circle = generate_circle(center, radius, 100)
# idxs, dists = find_nearest_neighbors(xs_circle, ys_circle, xs_noisy, ys_noisy)
# plot_projection(xs_circle, ys_circle, xs_noisy, ys_noisy, idxs)

