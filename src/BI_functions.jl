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


function llhood_circle(θs, model, xs_noisy, ys_noisy, σ)
    center = (θs[1], θs[2])
    radius = θs[3]
    # Compute radial deviations
    deltas = [sqrt((x - center[1])^2 + (y - center[2])^2) - radius for (x, y) in zip(xs_noisy, ys_noisy)]
    N = length(deltas)
    # Multivariate normal with mean zero and covariance σ^2 I
    mvn = MvNormal(zeros(N), Diagonal(fill(σ^2, N)))
    return logpdf(mvn, deltas)
end

function llhood_PD(θs, modelfn, neighboursfn,σn, σbias, ρij)
    # # scale = exp(θs[end])
    # # m     = modelfn(θs[1:end-1])
    
    # m = modelfn(θs)

    # if isempty(m)
    #     return -Inf
    # end
    
    rm, zm = modelfn(θs)

    if isempty(rm) || isempty(zm)
        return -Inf
    end

    m = hcat([rm, zm]...)  

    idxs, deltas    = neighboursfn(m)

    N               = length(deltas)
    # mvn             = MvNormal(zeros(N), Diagonal(fill(σn^2, N))) # no correlation

    # Σtry            = make_pd_cov(σn^2*ρij)
    # mvn             = MvNormal(zeros(N), Σtry) # including correlation

    mvn             = MvNormal(zeros(N), σn^2*ρij) # including correlation

    # mvn             = MvNormal(zeros(N), Diagonal(fill(σn^2, N)) + (scale*σbias)^2*ρij)
    return logpdf(mvn, vcat(deltas...))  

end

# make PD and symmetric by adding tiny diagonal jitter if necessary
function make_pd_cov(S::AbstractMatrix{<:Real}; eps0::Float64=1e-12, max_tries::Int=12)
    # symmetrize first
    Ssym = Symmetric((S + S') / 2)
    eps = eps0
    n = size(Ssym, 1)
    for i in 1:max_tries
        try
            # check via cholesky; raises if not PD
            cholesky(Ssym + eps * I(n); check=true)
            return Ssym + eps * I(n)
        catch e
            # increase jitter and retry
            eps *= 10.0
        end
    end
    error("make_pd_cov: unable to make matrix positive definite (tried up to eps=$(eps))")
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

###################################
### transformation to log space ###
###################################

# using mean, std
function mustdln(mu, sigma)

    muln = log(mu^2 / sqrt(mu^2 + sigma^2))
    stdln = sqrt(log(1 + sigma^2 / mu^2))
    return muln, stdln
end

# using mode, std
function mostdln(modeX::Float64, sigmaX::Float64)
    @assert modeX > 0 "modeX must be > 0"
    @assert sigmaX >= 0 "sigmaX must be >= 0"

    s2_over_m02 = (sigmaX^2) / (modeX^2)

    f(t) = t^4 - t^3 - s2_over_m02   # root for t > 1

    # bracket for t: must be > 1. Use (1+eps, big) but pick a practical upper bound
    a = 1.0000001
    b = 1e6   # very large upper bound; adjust if needed

    # Check sign change on interval; if not, try increasing b or fail
    if f(a)*f(b) > 0
        error("No sign change found for root search — inputs may be inconsistent.")
    end

    t = find_zero(f, (a,b), Bisection())  # robust bracketing solver
    σ2 = log(t)
    σ = sqrt(σ2)
    μ = log(modeX) + σ2

    # derived linear moments (for checking)
    meanX = exp(μ + 0.5*σ2)                 # or: modeX * exp(1.5*σ2)
    sigmaX_check = sqrt((exp(σ2)-1)*exp(2μ + σ2))

    skewX = (exp(σ2) + 2) * sqrt(exp(σ2) - 1)

    # return (mu=μ, sigma=σ, meanX=meanX, sigmaX_check=sigmaX_check, sigma2=σ2, t=t)
    return μ, σ, meanX, skewX
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


####################################
### llhood data normal direction ###
####################################



function bias_covariance_normal_local(X::AbstractMatrix, N::AbstractMatrix;
                                      σbias::Real=0.3, lbias::Real=1.0)
    n, d = size(X)
    @assert size(N) == (n, d) "Normals must be n×d (same number of points and dimensions)"

    # Normalize all normals
    Nnorm = similar(N)
    for i in 1:n
        ni = @view N[i, :]
        Nnorm[i, :] .= ni ./ sqrt(sum(abs2, ni))
    end

    Σ = Matrix{Float64}(undef, n, n)
    for i in 1:n
        Σ[i, i] = σbias^2
        xi = @view X[i, :]
        ni = @view Nnorm[i, :]
        for j in i+1:n
            xj = @view X[j, :]
            nj = @view Nnorm[j, :]
            nij = ni .+ nj
            denom = sqrt(sum(abs2, nij))
            nuse = denom > 1e-12 ? (nij ./ denom) : ni
            dij_n = abs(sum((xi .- xj) .* nuse))  # projection along averaged normal
            ρij = exp(-dij_n / lbias)
            val = σbias^2 * ρij
            Σ[i, j] = val
            Σ[j, i] = val
        end
    end
    return Σ
end



function bias_covariance(γ̇, σbias)
    # Correlation function based on distance between γ̇ values
    lbias = 1.0
    ρ(γ̇₁, γ̇₂) = exp(-abs(γ̇₁ - γ̇₂) / lbias)
    
    # Create covariance matrix with correlation structure
    n = length(γ̇)
    Σ = zeros(n, n)
    for i in 1:n
        for j in 1:n
            Σ[i, j] = σbias^2 * ρ(γ̇[i], γ̇[j])
        end
    end
    
    return Σ
end

function find_nearest_neighbors(m, xs_noisy, ys_noisy)
    xs_circle = m[:,1]
    ys_circle = m[:,2]
    perfect_points = hcat(xs_circle, ys_circle)'  # 2×N
    noisy_points = hcat(xs_noisy, ys_noisy)'      # 2×M
    tree = KDTree(perfect_points)
    idxs, dists = knn(tree, noisy_points, 1)      # 1 nearest neighbor per noisy point
    return getindex.(idxs, 1), dists  # Convert vector of vectors to flat vector
end


function find_nearest_neighbors2(m, xs_noisy, ys_noisy)
    xs_circle = m[:, 1]
    ys_circle = m[:, 2]
    perfect_points = hcat(xs_circle, ys_circle)'  # 2×N
    noisy_points   = hcat(xs_noisy, ys_noisy)'    # 2×M

    tree = KDTree(perfect_points)
    idxs, dists = knn(tree, noisy_points, 1)      # 1 nearest neighbor per noisy point
    idxs = getindex.(idxs, 1)                     # flatten indices

    # Compute direction vectors: noisy - model
    directions = [ (noisy_points[:, i] - perfect_points[:, idxs[i]]) for i in 1:length(idxs) ]

    # Optionally normalize to unit vectors
    unit_directions = [ d / norm(d) for d in directions ]

    return idxs, dists, directions, unit_directions
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


### compute ppd
function compute_ppd(Nsamples, modelfn, flatchain, rd, zd, Σ; norsamples::Int=1)
    
    QoI             = []
    QoIppd          = []
    Wos             = []
    μWos            = []
    flatchainlin    = exp.(flatchain)

    coords  = [(rd[j], zd[j]) for j in 1:size(rd)[1]]
    volume  = axisymmetric_volume_right(coords)

    for i = 1:Nsamples

        rm, zm              = modelfn(flatchain[:,i])
        Y                   = hcat([rm, zm]...)  
        idxs, deltas        = find_nearest_neighbors(Y, rd, zd)

        Ym = Y[idxs,:]
        rm = Ym[:,1]
        zm = Ym[:,2]

        normals = compute_normals(Ym, 10)
        # ptl = 0.9/128
        # σn = ptl / 2
        samples = sample_along_normals(Ym, normals, Σ, norsamples)
        ppdsam  = vcat(samples'...)

        wo      = wo_fn(volume, flatchainlin[1,end-i])

        push!(QoI, Ym)
        push!(QoIppd, ppdsam)
        push!(Wos, wo)
    end

    μγ      = mean(flatchainlin[1,:])
    μWo     = wo_fn(volume, μγ)

    return QoI, QoIppd, Wos, μWo
end

function QoIppd_plotting(Nsamples, samplesppd, QoIppd, flatchain)
    for i in 1:Nsamples
        samplesppd[i] = QoIppd[i] .- QoIppd[i][end,2]
    end


    meanθs = mean(flatchain, dims=2)
    rmmean, zmmean = modelfn(vec(meanθs))
    meanshape = hcat([rmmean, zmmean]...)

    alignedpoint = align_to_reference(meanshape, samplesppd)
    normalsref = compute_normals_2d(meanshape)
    lower_curve, upper_curve = credible_band(meanshape, normalsref, alignedpoint)

    x_mean = meanshape[:, 1]
    y_mean = meanshape[:, 2]

    x_upper = upper_curve[:, 1]
    y_upper = upper_curve[:, 2]

    x_lower = lower_curve[:, 1]
    y_lower = lower_curve[:, 2]

    return x_mean, y_mean, x_upper, y_upper, x_lower, y_lower

end