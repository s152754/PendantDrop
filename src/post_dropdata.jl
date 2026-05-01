using JLD2

# include("src/PDivp_model.jl")
# include("src/calculate_area_volume.jl")

function determine_arclength(rdfull::Vector{Float64}, zdfull::Vector{Float64})
    n = length(rdfull)
    ds = [sqrt((rdfull[i+1] - rdfull[i])^2 + (zdfull[i+1] - zdfull[i])^2) for i in 1:n-1]
    sdfull = vcat(0.0, cumsum(ds))  # cumulative arc length starting at 0
    return ds, sdfull
end

function trim_data(Nrep::Int, rdfull::Vector{Float64}, zdfull::Vector{Float64}, sdfull::Vector{Float64})
    # Target arc lengths evenly spaced
    s_target = range(first(sdfull), last(sdfull), length=Nrep)
    
    # Find closest indices for each target arc length
    idxs = [findmin(abs.(sdfull .- s))[2] for s in s_target]
    
    # Extract corresponding points
    rd = rdfull[idxs]
    zd = zdfull[idxs]
    sd = sdfull[idxs]
    
    return rd, zd, sd
end

# # exp. data
# Lp = load("data/droplet5/dropletcoords_rhs.jld2")
# rdf = Lp["rr_sorted"]
# zdf = Lp["Zr_sorted"]

# rdfull, zdfull, lastval = trim_trailing_repeats(rdf, zdf; keepmissing=false)
# R_cap = rd[end]
# ds, sdfull = determine_arclength(rdfull, zdfull)


# Nreps = [2, 4, 8, 16, 32, 64, 128]#, 256, 512]
# figdat = scatter(rdfull, zdfull, label="full", aspect_ratio=:equal)
# figrep = plot(xlabel="size data", ylabel="sdpart / sdfull", xaxis=:log)
# for i in eachindex(Nreps)
#     # i = 1
#     Nrep = Nreps[i]
    
#     rd, zd = trim_data(Nrep, rdfull, zdfull, sdfull)

#     fig = scatter(rdfull, zdfull, label="all data", aspect_ratio=:equal)
#     scatter!(rd, zd, label="partial data")
#     display(fig)

#     dspart, sdpart = determine_arclength(rd, zd)
#     scatter!(figrep, [(length(dspart)+1) / length(rdfull)], [(sdfull[end] - sdpart[end]) / sdfull[end]])

#     if i == 4
#         scatter!(figdat, rd, zd, label="partial")
#     end

# end

# display(figrep)



# display(figdat)



function trim_top_duplicates(x::Vector{Float64}, y::Vector{Float64}; tol=1e-12)
    n = length(x)
    
    # Start from top-left going down
    left_idx = 1
    while left_idx < n && abs(x[left_idx] - x[1]) < tol
        left_idx += 1
    end
    # Keep one point before change
    left_keep = max(left_idx - 1, 1)
    
    # Start from top-right going down (end of array)
    right_idx = n
    while right_idx > 1 && abs(x[right_idx] - x[end]) < tol
        right_idx -= 1
    end
    # Keep one point before change
    right_keep = min(right_idx + 1, n)
    
    # Build new arrays
    # x_new = vcat(x[1:left_keep], x[left_idx:right_keep], x[right_keep:n])
    # y_new = vcat(y[1:left_keep], y[left_idx:right_keep], y[right_keep:n])
    x_new = x[left_keep:right_keep]
    y_new = y[left_keep:right_keep  ]
    
    return x_new, y_new
end



# """
# Align each sample shape to the reference shape using Procrustes alignment.
# Inputs:
# - ref :: Matrix{Float64} (2 × N): reference shape points
# - samples :: Vector{Matrix{Float64}}: list of sample shapes (2 × N each)
# Returns:
# - aligned_samples :: Vector{Matrix{Float64}}
# """
# function align_to_reference(ref::Matrix{Float64}, samples::Vector{Matrix{Float64}})
#     aligned = Vector{Matrix{Float64}}(undef, length(samples))
#     for (i, s) in enumerate(samples)
#         # Center shapes
#         s_centered = s .- mean(s, dims=2)
#         ref_centered = ref .- mean(ref, dims=2)

#         # Compute optimal rotation via SVD
#         U, _, Vt = svd(s_centered * ref_centered')
#         R = U * Vt

#         # Apply rotation
#         aligned[i] = R * s_centered
#     end
#     return aligned
# end



# """
# Compute unit normals for a 2D curve.
# Input:
# - curve :: Matrix{Float64} (2 × N): points in order
# Returns:
# - normals :: Matrix{Float64} (2 × N): unit normals at each point
# """
# function compute_normals_2d(curve::Matrix{Float64})
#     N = size(curve, 2)
#     normals = zeros(2, N)
#     for i in 1:N
#         # Forward difference for tangent
#         j = i < N ? i + 1 : i - 1
#         tangent = curve[:, j] - curve[:, i]
#         tangent /= norm(tangent)
#         # Rotate tangent by 90° to get normal
#         normals[:, i] = [-tangent[2], tangent[1]]
#     end
#     return normals
# end



# """
# Compute upper and lower credible bands along normals.
# Inputs:
# - ref :: Matrix{Float64} (2 × N): reference shape
# - normals :: Matrix{Float64} (2 × N): normals at each point
# - aligned_samples :: Vector{Matrix{Float64}}: aligned sample shapes
# Returns:
# - lower_curve, upper_curve :: Matrix{Float64} (2 × N each)
# """
# function credible_band(ref::Matrix{Float64}, normals::Matrix{Float64}, aligned_samples::Vector{Matrix{Float64}})
#     N = size(ref, 2)
#     M = length(aligned_samples)
#     offsets = zeros(M, N)
#     for (m, sample) in enumerate(aligned_samples)
#         for i in 1:N
#             # Project displacement onto normal
#             disp = sample[:, i] - ref[:, i]
#             offsets[m, i] = dot(disp, normals[:, i])
#         end
#     end
#     # Compute quantiles for credible band (e.g., 95%)
#     lower = zeros(2, N)
#     upper = zeros(2, N)
#     for i in 1:N
#         q_low = quantile(offsets[:, i], 0.025)
#         q_high = quantile(offsets[:, i], 0.975)
#         lower[:, i] = ref[:, i] + q_low * normals[:, i]
#         upper[:, i] = ref[:, i] + q_high * normals[:, i]
#     end
#     return lower, upper

# end





# ### created filled droplet
# case = 5
# Ltot = load("data/droplet$(case)/dropletcoords.jld2")
# rtotfull = Ltot["rs_sorted"]
# ztotfull = Ltot["Zs_sorted"]



# xnew, ynew = reorder_droplet(rtotfull, ztotfull)
# xsnew, ysnew = trim_top_duplicates(xnew, ynew)

# fignew = plot(aspect_ratio=:equal) #xnew, ynew, 
# scatter!(fignew, xsnew, ysnew)
# display(fignew)

# x_closed = vcat(xsnew, xsnew[1])
# y_closed = vcat(ysnew, ysnew[1])

# figfill = plot(x_closed, y_closed, seriestype=:shape, color=:gray)
# plot!(figfill, aspect_ratio=:equal, grid=false, xaxis=false, yaxis=false, legend=false)
# display(figfill)
# savefig(figfill, "data/filled/droplet$(case).png")