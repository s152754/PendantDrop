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

#######################
### include scripts ###
#######################

include("src/BI_functions.jl")
include("src/PD_model.jl")
include("src/PDivp_model.jl")
include("src/io.jl")
include("src/utils.jl")
include("src/post_dropdata.jl")



# function sample_uniform(v::AbstractVector, N::Int)
#     len = length(v)
#     if N > len
#         error("N cannot be greater than the length of the vector")
#     end
#     # Compute indices uniformly spaced
#     idx = round.(Int, range(1, len; length=N))
#     return v[idx]
# end

case = 3

###################
### import data ###
###################

# real data
Lp = load("data/droplet$(case)/dropletcoords_rhs_out.jld2")
rdf = Lp["rr_sorted"]
zdf = Lp["Zr_sorted"]

rdfull, zdfull, lastval = trim_trailing_repeats(rdf, zdf; keepmissing=false)

Nrep = 16 #Int(round(0.05*length(rdfull)))
ds, sdfull = determine_arclength(rdfull, zdfull)
rd, zd, sd = trim_data(Nrep, rdfull, zdfull, sdfull)
# rd = rdfull
# zd = zdfull
# sd = sdfull

# rdfull, zdfull, lastval = trim_trailing_repeats(rdf, zdf; keepmissing=false)
# Nrep = Int(round(0.1*length(rdfull)))
# rd = sample_uniform(rdfull,Nrep)
# zd = sample_uniform(zdfull, Nrep)

figdatasam = plot(xlabel="r [mm]", ylabel="z [mm]", aspect_ratio=:equal, legend=false)
# scatter!(figdatasam, rdfull, zdfull, label="full data", color=:black)
scatter!(figdatasam, rd, zd, label="partial data", color=:red)
display(figdatasam)

Rneedle = rd[end]  # capillary radius in mm

#######################
### determine noise ###
#######################

ptl = 0.9/139
# wb = 7 # width boundary
# pxbound = (wb * ptl) # assume broad normal distribution
# wN = 6 # width of normal distribution
# σbias = pxbound / wN

# lϵs = [1e-12, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1, 1.]

# figcor = plot(xlabel="lag h", ylabel="ρ(h)")
# for i in eachindex(lϵs)
#     ρij = autocornoise(sd, lϵs[i])
#     plot!(figcor, 1:length(rd), ρij[1,:], label="l_ϵ=$(lϵs[i])", linestyl=:solid, marker=:circle)
# end
# display(figcor)

lϵ = 1e-12
ρij = autocornoise(sd, lϵ)
σn = 3*ptl  # noise is 0.5 pixel (ptl/2) #0.05*rd[end] # measurement noise in mm. dus 3 pixels total normal dist

##############################
### import prior and model ###
##############################
sigma = 1000.
θgum = 1.
σpr = 0.25
# Rcaps = [0.40, 0.50, 0.60, 0.65, 0.70]
Rcaps = [0.72, 1.1, 1.2, 1.3, 1.4]
Rcap = Rcaps[case]
# Rcaps = [2.28, 2.33, 2.47, 2.59, 2.67]
# Rcap = 1/Rcaps[case]

prior           = OrderedDict{String, Distribution}()
prior["γ"]      = Normal(mostdln(72.0, σpr*72.0)[1], mostdln(72.0, σpr*72.0)[2])
prior["Rap"]    = Normal(mostdln(Rcap, σpr*Rcap)[1], mostdln(Rcap, σpr*Rcap)[2])
# # prior["scale"]  = -Gumbel(-log(sigma), θgum)

# prior         = OrderedDict{String, Distribution}()
# prior["γ"]      = Normal(mostdln(72.0, 0.01*72.0)[1], mostdln(72.0, 0.01*72.0)[2])
# prior["Rap"]    = Normal(mostdln(Rcaps[5], σpr*Rcaps[5])[1], mostdln(Rcaps[5], σpr*Rcaps[5])[2])
# prior["scale"]  = -Gumbel(-log(sigma), θgum)

# model

modelfn         = θs -> pendant_drop_shape_RK4(θs, Rneedle, 5000, r_tol=1e-12, h_tol=1e-12, hits_max=2)# returns [r,z]
neighboursfn    = ms -> find_nearest_neighbors(ms, rd, zd)
llhoodfn        = θs -> llhood_PD(θs, modelfn, neighboursfn, σn, σbias, ρij)
lpriorfn        = θs -> lprior(θs, values(prior))
lprobfn         = θs -> lprob(θs, lpriorfn, llhoodfn)

#####################
### AISM sampling ###
#####################
numdims                 = length(prior)
numwalkers              = 4*numdims
numsamples_perwalker    = 1500 # per walker
burnin                  = 800 # per walker
thinning                = 1
a                       = 2.0
global burncount        = numwalkers*(burnin + 1)


# priorμ                  = log.([72., Rcaps[1], 1.1])
priorμ                  = log.([72., Rcap])
priorσ                  = 0.001
ps                      = dictsamplemean(priorμ, abs.(priorμ)*priorσ, numwalkers) 
# ps                          = Matrix(dictsample(prior, numwalkers)') # sample from prior
chain, llhoodvals           = AISMburn(lprobfn, numwalkers, ps, burnin,1) # burnin time
chain, llhoodvals           = AISM(lprobfn, numwalkers, chain[:, :, end], numsamples_perwalker, thinning, a) # sampling after burnin

flatchain, flatllhoodvals   = AISMflatten(chain, llhoodvals)
Nsamples = numsamples_perwalker*numwalkers 

jldsave("data/droplet$(case)/chain_out_final.jld2"; chain, llhoodvals, flatchain, flatllhoodvals)

###################
### chain plots ###
###################

flatchainlin = exp.(flatchain)

# post vs prior plot
n_cols = ceil(Int, sqrt(numdims))
n_rows = ceil(Int, numdims / n_cols)

# Maak subplots dynamisch
p = plot(layout=(n_rows, n_cols), legend=false)
labels = collect(keys(prior))

for i in eachindex(labels)
        xtick_positions = range(minimum(flatchainlin[i,:]), stop=maximum(flatchainlin[i,:]), length=3)
        xtick_labels = [string(round(x, sigdigits=3)) for x in xtick_positions]
        histogram!(p, flatchainlin[i,:], xlabel=labels[i], subplot=i, normalize=:pdf, color=:black, xticks = (xtick_positions, xtick_labels), label="")
        xp = range(minimum(flatchain[i,:]), stop=maximum(flatchain[i,:]), length=200)
        xplin = exp.(xp)
        J = exp.(xp)
        plot!(p, xplin, pdf.(collect(values(prior))[i], xp)./J, subplot=i,color=:red, lw=2, label="")
end
display(p)      

# trace plot
pt = plot(layout=(n_rows, n_cols), legend=false)
for i in eachindex(labels)
    plot!(pt, chain[i,:,:]', xlabel=labels[i], subplot=i)
end
display(pt)

nplot = 100
figcomp = plot(aspect_ratio=:equal)
QoI, QoIppd      = compute_ppd(nplot, modelfn, flatchain, rd, zd, σn)
for i = 1:nplot
# i = 90
    # Y = modelfn(flatchain[:,end - i])
    # r = Y[:, 1]
    # z = Y[:, 2]

    r, z = modelfn(flatchain[:,end - i])

    plot!(figcomp, r, z .- z[end]; alpha=0.05, color=:red, label="")#, ylims=(-2.0,0)
    
    # plot!(figcomp, QoI[i][:,1], QoI[i][:,2], alpha=0.5, color=:blue, label="")
end
scatter!(figcomp, rd, zd, label="exp", color=:black, aspect_ratio=:equal, markersize=1)
display(figcomp)


println("mean = $(mean(flatchainlin, dims=2)), mode = $(mode(flatchainlin[1,:])), std = $(std(flatchainlin, dims=2)), CV = $(std(flatchainlin, dims=2) ./ mean(flatchainlin, dims=2) * 100 )")








## compute real ppd for figures
using Statistics
include("src/calculate_area_volume.jl")
include("src/BI_functions.jl")
ptl         = 0.9/139 # 128
σn          = 3*ptl
Σ           =  σn
Nsamples    = 2000
ncase = 5



L = load("data/droplet$(ncase)/chain_out_final.jld2")

flatchain = L["flatchain"]

# exp. data
Lp = load("data/droplet$(ncase)/dropletcoords_rhs_out.jld2")
rdf = Lp["rr_sorted"]
zdf = Lp["Zr_sorted"]

rdfull, zdfull, lastval = trim_trailing_repeats(rdf, zdf; keepmissing=false)

Nrep        = 16
ds, sdfull  = determine_arclength(rdfull, zdfull)
rd, zd, sd  = trim_data(Nrep, rdfull, zdfull, sdfull)
Rneedle     = rd[end]

modelfn         = θs -> PD_fn(θs, Rneedle, n_points=5000, r_tol=1e-12, h_tol=1e-12, hits_max=2)# returns [r,z]
QoI, QoIppd, Wo, μWo = compute_ppd(Nsamples, modelfn, flatchain, rd, zd, Σ)

jldsave("output/QoI_droplet$(ncase).jld2"; QoI, QoIppd, Wo, μWo)

L = load("output/QoI_droplet$(ncase).jld2")["μWo"]
println("$(L)")
##






# #########
# ## ppd ##
# #########

# QoI = []
# QoIppd = []
# for i = 1:100
#     # i = 10
#     Y = modelfn(flatchain[:,i])

#     idxs, deltas = find_nearest_neighbors(Y, rd, zd)

#     Ym = Y[idxs,:]
#     rm = Ym[:,1]
#     zm = Ym[:,2]

#     normals = compute_normals(Ym, 10)
#     ptl = 0.9/128
#     σn = ptl / 2
#     numsamples = 1
#     samples = sample_along_normals(Ym, normals, σn, numsamples)
#     ppdsam = vcat(samples'...)
#     push!(QoI, Ym)
#     push!(QoIppd, ppdsam)
# end
# # QoI = hcat(QoItest...)
# # QoIppd = hcat(QoIppdtest...)

# # μQoI      = mean(QoI,dims=2)
# # σQoI      = std(QoI,dims=2)
# # μQoIppd   = mean(QoIppd, dims=2)
# # σQoIppd   = std(QoIppd, dims=2)


# figdatas = scatter(rd, zd; label="exp data", color=:blue, markersize=2, aspect_ratio=:equal, xlabel="r (mm)", ylabel="z (mm)")
# scatter!(rdf, zdf .- zdf[end]; label="exp data (full)", color=:gray, markersize=1, alpha=0.3)
# display(figdatas)
