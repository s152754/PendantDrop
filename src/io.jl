function trim_trailing_repeats(r::Vector{T}, z::Vector{S}; keepmissing::Bool = false) where {T,S}
    length(r) == length(z) || throw(ArgumentError("r and z must have the same length"))

    n = length(r)
    if n == 0
        return r, z, nothing
    end

    # Find the last non-missing value in r
    lastidx = nothing
    for i in n:-1:1
        if keepmissing ? true : r[i] !== missing
            lastidx = i
            if r[i] === missing && !keepmissing
                continue
            end
            break
        end
    end

    if lastidx === nothing || r[lastidx] === missing
        return r, z, nothing
    end

    lastval = r[lastidx]

    # Walk upwards from the end while r[i] == lastval
    cut = n
    while cut >= 1
        if r[cut] === lastval
            cut -= 1
        else
            break
        end
    end

    # Keep one repeated value (cut+1)
    r2 = r[1:cut+1]
    z2s = z[1:cut+1]

    # Shift z so last value becomes 0
    z2 = z2s .- z2s[end]

    return r2, z2, lastval
end

function align_to_reference(reference::Matrix{Float64}, samples::Vector{Matrix{Float64}})
    N = size(reference, 1)
    aligned_points = [zeros(N, 2) for _ in 1:length(samples)]
    trees = [KDTree(s') for s in samples]

    for i in 1:length(samples)
        tree = trees[i]
        for j in 1:N
            idx, _ = knn(tree, reference[j, :], 1, true)
            aligned_points[i][j, :] = samples[i][idx[1], :]
        end
    end

    return aligned_points
end



function compute_normals_2d(curve::Matrix{Float64})
    N = size(curve, 1)
    normals = zeros(N, 2)
    for i in 1:N
        if i == 1
            t = curve[i+1, :] .- curve[i, :]
        elseif i == N
            t = curve[i, :] .- curve[i-1, :]
        else
            t = curve[i+1, :] .- curve[i-1, :]
        end
        n = [-t[2], t[1]]
        normals[i, :] = n / norm(n)
    end
    return normals
end



function credible_band(reference::Matrix{Float64}, normals::Matrix{Float64}, aligned_points::Vector{Matrix{Float64}}; alpha=0.05)
    N = size(reference, 1)
    lower_curve = zeros(N, 2)
    upper_curve = zeros(N, 2)

    for j in 1:N
        n = normals[j, :]
        dists = [dot(aligned_points[i][j, :] .- reference[j, :], n) for i in 1:length(aligned_points)]
        lower_q = quantile(dists, alpha/2)
        upper_q = quantile(dists, 1 - alpha/2)
        lower_curve[j, :] = reference[j, :] .+ lower_q .* n
        upper_curve[j, :] = reference[j, :] .+ upper_q .* n
    end

    return lower_curve, upper_curve
end

"""
Compute normals for a point cloud using local PCA.
points: Nx2 matrix of droplet points
k: number of neighbors for PCA
"""
function compute_normals(points::Matrix{Float64}, k::Int=10)
    N = size(points, 1)
    normals = zeros(N, 2)
    
    # Build KDTree (points' columns are dimensions)
    tree = KDTree(points')
    
    for i in 1:N
        # Find k nearest neighbors (includes point i)
        idxs, _ = knn(tree, points[i, :], k, true)
        neighbors = Array(points[idxs, :])  # ensure a concrete matrix

        # skip if not enough neighbors to form covariance
        if size(neighbors, 1) < 2
            normals[i, :] .= [0.0, 0.0]
            continue
        end

        # fully-qualified covariance to avoid ambiguity
        covmat = Statistics.cov(neighbors)  # returns 2x2 if neighbors is k×2

        # eigen decomposition
        eig = eigen(Symmetric(covmat))      # ensure symmetric
        # eigenvalues in eig.values, eigenvectors in eig.vectors (columns)
        normal = eig.vectors[:, argmin(eig.values)]

        # normalize (just in case) and store as row
        normals[i, :] = (normal / norm(normal))'
    end
    
    return normals
end


"""
Generate samples along normal direction for each point.
points: Nx2 matrix
normals: Nx2 matrix
sigma: standard deviation along normal
num_samples: number of samples per point
"""
function sample_along_normals(points::Matrix{Float64}, normals::Matrix{Float64}, sigma::Float64, num_samples::Int=1)
    N = size(points, 1)
    samples = Vector{Matrix{Float64}}(undef, N)  # each entry: 2 x num_samples

    for i in 1:N
        p = reshape(points[i, :], 2, 1)                 # 2×1
        ncol = reshape(normals[i, :], 2, 1)
        if norm(ncol) == 0.0
            # if normal undefined, produce repeated original point(s)
            samples[i] = repeat(p, 1, num_samples)
            continue
        end
        ncol ./= norm(ncol)                             # unit normal (2×1)
        δs = randn(num_samples) .* sigma                # 1×num_samples
        samples[i] = p .+ ncol * δs'                    # 2×num_samples
    end

    return samples
end


function get_data(dir_data, filename; Nreps=16)
    isdir(dir_data) || error("Directory not found: $dir_data")

    file = joinpath(dir_data, filename)
    @load file rs_sorted Zs_sorted

    rptots, zptots = reorder_droplet(rs_sorted, Zs_sorted)
    rptotf, zptotf = trim_top_duplicates(rptots, zptots)

    rd, zd, sd = trim_data(Nreps, rptotf, zptotf)

    return rd, zd
end