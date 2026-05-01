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

#######################
### include scripts ###
#######################

include("../src/BI_functions.jl")

####################
### circle model ###
####################

function generate_circle(params, n_points)
    x_c = params[1]
    y_c = params[2]
    radius = params[3]

    xs = [x_c + radius * cos(θ) for θ in range(0, 2π, length=n_points)]
    ys = [y_c + radius * sin(θ) for θ in range(0, 2π, length=n_points)]
    return hcat(xs, ys)
end

#########################
### noisy data points ###
#########################

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

##############################
### find nearest neighbors ###
##############################
function find_nearest_neighbors(m, xs_noisy, ys_noisy)
    xs_circle = m[:,1]
    ys_circle = m[:,2]    
    perfect_points = hcat(xs_circle, ys_circle)'  # 2×N
    noisy_points = hcat(xs_noisy, ys_noisy)'      # 2×M
    tree = KDTree(perfect_points)
    idxs, dists = knn(tree, noisy_points, 1)      # 1 nearest neighbor per noisy point
    return getindex.(idxs, 1), dists  # Convert vector of vectors to flat vector
end

#############
### prior ###
#############
prior           = OrderedDict{String, Distribution}()
prior["x_c"]    = Normal(1.,0.1*1.0)
prior["y_c"]    = Normal(-0.5, 0.1*0.5)
prior["r"]      = Normal(5., 0.1*5.0)


####################
### start script ###
####################
xcm = 0.0
ycm = 0.0
center = (xcm, ycm)
radius = 5.0
σn = 0.1*radius
xs_circle, ys_circle = generate_circle([xcm, ycm, radius], 100)
xs_noisy, ys_noisy = generate_noisy_points(center, radius, 50; noise_std=σn)


# bi functions
npoints         = 100
modelfn         = θs -> generate_circle(θs, npoints)  # returns (center, radius)
# llhoodfn        = θs -> llhood_circle(θs, modelfn, xs_noisy, ys_noisy, σn)
neighboursfn    = ms -> find_nearest_neighbors(ms, xs_noisy, ys_noisy)
llhoodfn        = θs -> llhood_PD(θs, modelfn, neighboursfn, σn)
lpriorfn        = θs -> lprior(θs, values(prior))
lprobfn         = θs -> lprob(θs, lpriorfn, llhoodfn)


# aism sampling

numdims                 = length(prior)
numwalkers              = 4*numdims # 32*numdims
numsamples_perwalker    = 1000 #1250
burnin                  = 1 # per walker
thinning                = 1
a                       = 2.0
global burncount        = numwalkers*(burnin + 1)

ps                          = Matrix(dictsample(prior, numwalkers)') # sample from prior
chain, llhoodvals           = AISMburn(lprobfn, numwalkers, ps, burnin,1) # burnin time
# println("relative excluded number of paramsets [0,1] = $(NOTpossible/(size(chain)[3]*numwalkers))")
chain, llhoodvals           = AISM(lprobfn, numwalkers, chain[:, :, end], numsamples_perwalker, thinning, a) # sampling after burnin

flatchain, flatllhoodvals   = AISMflatten(chain, llhoodvals)
Nsamples = numsamples_perwalker*numwalkers 

# QoI
QoItest = []
for i = 1:3#size(flatchain)[2]
    qoitest1, qoitest2 = modelfn(flatchain[:,i])
    qoitest = [qoitest1, qoitest2]
    println("qoi = $(qoitest[end, end])")
    push!(QoItest, qoitest)
end

QoI = hcat(QoItest...)

# post processing
# figcircle = scatter(xs_noisy, ys_noisy, label="exp", color=:black, aspect_ratio=:equal)
# scatter!(figcircle, xs_circle, ys_circle, label="model")
# # plot!(figcircle, , QoI[:,end-nplot:end], alpha=0.05, ylims=(0,400), color=:red, label="")
# display(figcircle)

# post vs prior plot
n_cols = ceil(Int, sqrt(numdims))
n_rows = ceil(Int, numdims / n_cols)

# Maak subplots dynamisch
p = plot(layout=(n_rows, n_cols), legend=false)
labels = collect(keys(prior))

# Fill subplots in a loop
for i in eachindex(labels)
    histogram!(p, flatchain[i,:], xlabel=labels[i], subplot=i, normalize=:pdf)
    xp = range(minimum(flatchain[i,:]), stop=maximum(flatchain[i,:]), length=200)
    plot!(p, xp, pdf.(collect(values(prior))[i], xp), subplot=i,color=:red, lw=2)
end

# Display plot
display(p)        

# trace plot
pt = plot(layout=(n_rows, n_cols), legend=false)
for i in eachindex(labels)
    plot!(pt, chain[i,:,:]', xlabel=labels[i], subplot=i)
end
display(pt)

