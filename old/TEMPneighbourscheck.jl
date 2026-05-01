using Plots
using NearestNeighbors
using Distributions, LinearAlgebra
using JLD2
using OrderedCollections

include("src/BI_functions.jl")
include("src/PDivp_model.jl")
include("src/io.jl")
include("src/post_dropdata.jl")

# real data
Lp = load("data/droplet5/dropletcoords_rhs.jld2")
rdf = Lp["rr_sorted"]
zdf = Lp["Zr_sorted"]

rdfull, zdfull, lastval = trim_trailing_repeats(rdf, zdf; keepmissing=false)

Nrep = 16
ds, sdfull = determine_arclength(rdfull, zdfull)
rd, zd = trim_data(Nrep, rdfull, zdfull, sdfull)

R_cap = rd[end]

modelfn         = θs -> pendant_drop_shape_RK4(θs, R_cap, 10000, r_tol=1e-12, h_tol=1e-12, hits_max=2)# returns [r,z]


rm, zm = modelfn(log.([72.0, 0.7]))  # example parameters: sigma, R
m = hcat([rm, zm]...)  
idxs, dists = find_nearest_neighbors(m, rd, zd)

figproj= plot(xlabel="r [mm]", ylabel="z [mm]", aspect_ratio=:equal, legend=:topleft)
# plot!(figproj, ylims=(-4.0,-3.0), xlims=(0,1.0))
plot!(figproj, rm, zm, label="model", color=:blue, aspect_ratio=:equal)
plot!(figproj, rd, zd, label="exp. data", color=:red)
# Highlight closest neighbors
xs_closest = rm[idxs]
ys_closest = zm[idxs]
for i in eachindex(rd)
    plot!(figproj, [rd[i], xs_closest[i]], [zd[i], ys_closest[i]], color=:gray, alpha=0.5, label="")
end
display(figproj)