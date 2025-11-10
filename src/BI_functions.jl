function lprior(θs, dists)
    logPrior = sum(logpdf(dist, θ) for (θ, dist) in zip(θs, dists))

    if logPrior == -Inf
        return -Inf
    end

    return logPrior
end

function llhood(θs, model, μdata, σdata)
    m = model(θs)

    if isempty(m)
        return -Inf
    end

	return logpdf(MvNormal(μdata, diagm(σdata.^2)), m)
end

function llhoodinfnoise(θs, model, μdata, σdata, ρ)
    scale = exp(θs[end])
    # scale = θs[end]
    m = model(θs[1:end-1])

    if isempty(m)
        return -Inf
    end

    σbias = scale*mean(σdata)

	return logpdf(MvNormal(μdata, diagm((σdata).^2)+σbias^2*ρ), m)
end

function autocornoise(tdata, lϵ)
    ρ = zeros(length(tdata), length(tdata))
    for i in eachindex(tdata)
        for j in eachindex(tdata)
            ρ[i,j] = exp(-abs.(tdata[i] - tdata[j])/lϵ)
        end
    end

    return ρ
end

function lprob(θs, lpriorfn, llhoodfn)    
    if isinf(lpriorfn(θs))
        return lpriorfn(θs)
    else
        return lpriorfn(θs) + llhoodfn(θs)
    end    
end

####################
### AISM package ###
####################

function AISMburn(llhood, numwalkers, x₀, burnin, thinning)
    chain, chainlnp = AffineInvariantMCMC.sample(llhood, numwalkers, x₀, burnin, thinning)
    return chain, chainlnp
end

function AISM(llhood, numwalkers, x₀, numsamples_perwalker, thinning, a)
    chain, chainlnp = AffineInvariantMCMC.sample(llhood, numwalkers, x₀, numsamples_perwalker, thinning, a)
    return chain, chainlnp
end

function AISMflatten(chain, chainlnp)
    flatchain, flatchainlnp = AffineInvariantMCMC.flattenmcmcarray(chain, chainlnp)
    return flatchain, flatchainlnp
end