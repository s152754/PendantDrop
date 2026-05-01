using Distributions, LinearAlgebra
using Statistics
using Plots
using JLD2
# using Turing
using Random
using AffineInvariantMCMC
# using SpecialFunctions
using OrderedCollections
using LaTeXStrings
using CSV
using NearestNeighbors
using Roots
using StatsBase

include("src/BI_functions.jl")
# include("src/PD_model.jl")
include("src/RK4model.jl")
include("src/io.jl")
include("src/post_dropdata.jl")

case = 5

# real data
Lp = load("data/droplet$(case)/dropletcoords_rhs_out.jld2")
rdf = Lp["rr_sorted"]
zdf = Lp["Zr_sorted"]

rdfull, zdfull, lastval = trim_trailing_repeats(rdf, zdf; keepmissing=false)


ds, sdfull = determine_arclength(rdfull, zdfull)
# Nrep = 100
# rd, zd, sd = trim_data(Nrep, rdfull, zdfull, sdfull)
lϵ = 1.



figcor = plot(xlabel="n data points", ylabel="correlation")
lϵs = range(0.1, 10.0, length=10)
for i = 1:10
    lϵ = lϵs[i]
    ρ = autocornoise(sdfull, lϵ)
    plot!(figcor, ρ[1,:], label="lϵ = $(lϵ))")
end

display(figcor)


# visualize droplet every 50th point
figdrop = plot(rdfull, zdfull, label="full", aspect_ratio=:equal)
plot!(figdrop, xlabel="r [mm]", ylabel="z [mm]")
scatter!(figdrop, rdfull[1:50:end], zdfull[1:50:end], color=:red, markersize=6, label="Every 50th point")
title!(figdrop, "every 50th point (from bottom) highlighted")
display(figdrop)