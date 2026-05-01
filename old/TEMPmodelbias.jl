using JLD2, Plots, OrderedCollections, Distributions, NearestNeighbors, LinearAlgebra, Random

include("src/io.jl")
include("src/BI_functions.jl")
include("src/PDivp_model.jl")


function check_bias_alignment(data::AbstractMatrix, model::AbstractMatrix, normals::AbstractMatrix)
    n, d = size(data)
    @assert size(model) == (n, d)
    @assert size(normals) == (n, d)

    ratios = zeros(n)
    angles = zeros(n)

    for i in 1:n
        # create explicit views then operate
        vi = @view(data[i, :])
        vj = @view(model[i, :])
        b = vi .- vj

        nvec = @view(normals[i, :])
        nunit = nvec ./ sqrt(sum(abs2, nvec))

        proj = sum(b .* nunit)
        norm_b = sqrt(sum(abs2, b))

        if norm_b == 0.0
            ratios[i] = 0.0
            angles[i] = 0.0
            continue
        end

        ratios[i] = abs(proj) / norm_b
        angles[i] = acos(clamp(proj / norm_b, -1.0, 1.0)) * (180 / π)  # degrees
    end

    return ratios, angles
end


# real data
Lp = load("data/droplet5/dropletcoords_rhs.jld2")
rdf = Lp["rr_sorted"]
zdf = Lp["Zr_sorted"]

rd, zd, lastval = trim_trailing_repeats(rdf, zdf; keepmissing=false)
R_cap = rd[end]

modelfn         = θs -> pendant_drop_shape_RK4(θs, R_cap, 10000, r_tol=1e-12, h_tol=1e-12, hits_max=2)# returns [r,z]


m = modelfn(log.([72.0, 0.7]))  # example parameters: sigma, R
rm = m[:,1]
zm = m[:,2]
idxs, dists, dirs, unit_dirs = find_nearest_neighbors2(m, rd, zd)
r = rm[idxs]
z = zm[idxs]
mm = hcat(r,z)
ptl = 0.9/128
σbias = ptl*3 # 3 pixels





# ------------------------------------------------------------------------------
# Helper: normalize normals provided as Vector{Vector{Float64}} → n×2 matrix
# ------------------------------------------------------------------------------
function normalize_normals(normal_vectors::Vector{Vector{Float64}})
    n = length(normal_vectors)
    N = Matrix{Float64}(undef, n, 2)
    for i in 1:n
        ni = normal_vectors[i]
        @assert length(ni) == 2 "Each normal must be length-2 [nx, ny]"
        norm_ni = sqrt(sum(abs2, ni))
        @assert norm_ni > 0 "Normal vector at index $i must be non-zero"
        N[i, 1] = ni[1] / norm_ni
        N[i, 2] = ni[2] / norm_ni
    end
    return N
end

# ------------------------------------------------------------------------------
# Build covariance Σ with correlation along local normals using lbias
# dn_ij = | (x_i - x_j) ⋅ n_ij |, with n_ij the average unit normal between i, j
# Kernel: :exp (Matérn ν=1/2-like) or :se (squared-exponential)
# ------------------------------------------------------------------------------
function covariance_along_normals(X::AbstractMatrix, N::AbstractMatrix;
                                  σbias::Real=0.3, lbias::Real=1.0, kernel::Symbol=:exp)
    n, d = size(X)
    @assert d == 2 "This implementation assumes 2D (x,y)"
    @assert size(N) == (n, d) "Normals must be n×2"

    ρ(dn) = kernel === :exp  ? exp(-abs(dn) / lbias) :
            kernel === :se   ? exp(-(dn^2) / (2 * lbias^2)) :
            throw(ArgumentError("kernel must be :exp or :se"))

    Σ = Matrix{Float64}(undef, n, n)
    # Fill symmetric matrix efficiently
    for i in 1:n
        Σ[i, i] = σbias^2
        xi = @view X[i, :]
        ni = @view N[i, :]
        for j in i+1:n
            xj = @view X[j, :]
            nj = @view N[j, :]

            nij = ni .+ nj
            denom = sqrt(sum(abs2, nij))
            # If normals nearly opposite, average can be tiny; fallback to ni
            nuse = denom > 1e-12 ? (nij ./ denom) : ni

            dn = abs(sum((xi .- xj) .* nuse))  # normal-projected separation
            val = σbias^2 * ρ(dn)
            Σ[i, j] = val
            Σ[j, i] = val
        end
    end
    return Σ
end

# ------------------------------------------------------------------------------
# Sample a correlated Gaussian vector with covariance Σ via Cholesky
# Add a tiny diagonal jitter to ensure numerical stability
# ------------------------------------------------------------------------------
function sample_correlated(Σ::AbstractMatrix; jitter::Real=1e-10, rng::AbstractRNG=Random.GLOBAL_RNG)
    n = size(Σ, 1)
    Σj = copy(Σ)
    @inbounds for i in 1:n
        Σj[i, i] += jitter
    end
    F = cholesky(Σj; check=false).L      # lower-triangular
    z = randn(rng, n)
    return F * z                          # b ~ N(0, Σj)
end

# ------------------------------------------------------------------------------
# Main: use lbias to create correlated bias magnitudes, then apply along normals
# ------------------------------------------------------------------------------
function model_plus_correlated_bias(model_points::Matrix{Float64},
                                    normal_vectors::Vector{Vector{Float64}};
                                    σbias::Real=0.3, lbias::Real=1.0,
                                    kernel::Symbol=:exp,
                                    rng::AbstractRNG = Random.MersenneTwister(1234))
    n, d = size(model_points)
    @assert d == 2
    @assert length(normal_vectors) == n

    # Normalize normals → n×2
    N = normalize_normals(normal_vectors)

    # Build Σ using positions (model_points) and local normals, with lbias
    Σ = covariance_along_normals(model_points, N; σbias=σbias, lbias=lbias, kernel=kernel)

    # Sample correlated magnitudes
    b = sample_correlated(Σ; rng=rng)   # length n vector

    # Apply along normals
    biased_points = similar(model_points)
    for i in 1:n
        biased_points[i, :] = model_points[i, :] .+ b[i] .* N[i, :]
    end
    return biased_points, b, Σ
end

# ------------------------------------------------------------------------------
# === Usage with your existing variables ===
# data_points :: Matrix{Float64} (n×2)
# model_points :: Matrix{Float64} (n×2)
# normal_vectors :: Vector{Vector{Float64}} (length n, each [nx, ny])
# ------------------------------------------------------------------------------
lbias = 1.0    # <<< change this to see the effect (e.g., 0.05, 0.5, 2.0)

biased_points, b, Σ = model_plus_correlated_bias(model_points, normal_vectors;
                                                 σbias=σbias, lbias=lbias, kernel=:exp)

# Print first 10 coordinates (as you requested earlier)
println("First 10 Experimental Points:\n", data_points[1:10, :])
println("\nFirst 10 Model Points:\n", model_points[1:10, :])
println("\nFirst 10 Biased Points (correlated by lbias):\n", biased_points[1:10, :])

# ------------------------------------------------------------------------------
# Plot: points + lines (no arrows)
# - Gray line: model → biased
# - Light gray line: model → experimental
# ------------------------------------------------------------------------------
plt = scatter(data_points[:, 1], data_points[:, 2], label="Experimental", color=:blue, markerstrokecolor=:blue)
scatter!(plt, model_points[:, 1], model_points[:, 2], label="Model", color=:green, markerstrokecolor=:green)
scatter!(plt, biased_points[:, 1], biased_points[:, 2], label="Model + Correlated Bias", color=:red, markerstrokecolor=:red)
plot!(plt, aspect_ratio=:equal, ylims=(-4, -3.25), xlims=(0, 1))
# Draw lines for first 10 points for clarity (increase if you want)
for i in 1:10
    msx, msy = model_points[i, :]
    bsx, bsy = biased_points[i, :]
    plot!(plt, [msx, bsx], [msy, bsy], color=:gray, linewidth=2, label="")
end
for i in 1:10
    msx, msy = model_points[i, :]
    dsx, dsy = data_points[i, :]
    plot!(plt, [msx, dsx], [msy, dsy], color=:lightgray, linewidth=2, label="")
end

title!("Effect of Autocorrelation Length lbias = $(lbias)")
xlabel!("r [mm]"); ylabel!("z [mm]")
display(plt)


# Σ = bias_covariance_normal_local(mm, vcat(unit_dirs'...), σbias=σbias, lbias=0.01)
# ratios, angles = check_bias_alignment(hcat(rd, zd), hcat(r, z), vcat(unit_dirs'...))

# println("Average alignment ratio: ", mean(ratios))
# println("Average angle (deg): ", mean(angles))

# # Combine into matrices
# data_points = hcat(rd, zd)
# model_points = hcat(r, z)
# normal_vectors = unit_dirs

# # -----------------------------
# # Generate biased points: model + bias along normal
# # -----------------------------
# Random.seed!(1234)  # reproducibility

# biased_points = similar(model_points)
# for i in 1:size(model_points, 1)
#     nvec = normal_vectors[i]
#     nunit = nvec ./ norm(nvec)  # normalize
#     bias_magnitude = randn() * σbias
#     biased_points[i, :] = model_points[i, :] .+ bias_magnitude .* nunit
# end

# # -----------------------------
# # Plot everything
# # -----------------------------
# figbias = scatter(data_points[:, 1], data_points[:, 2], label="Experimental", color=:blue)
# plot!(figbias, title="influence model bias", xlabel="r [mm]", ylabel="z [mm]", aspect_ratio=:equal, ylims=(-4, -3), xlims=(0, 1))
# scatter!(figbias, model_points[:, 1], model_points[:, 2], label="Model", color=:green)
# scatter!(figbias, biased_points[:, 1], biased_points[:, 2], label="Model + Bias", color=:red)


# # Draw light gray lines between model and experimental points
# for i in 1:10
#     x_start, y_start = model_points[i, :]
#     x_end, y_end = data_points[i, :]
#     plot!(figbias, [x_start, x_end], [y_start, y_end], color=:lightgray, linewidth=2, label="")
# end

# # Draw gray lines between model and biased points
# for i in 1:10
#     x_start, y_start = model_points[i, :]
#     x_end, y_end = biased_points[i, :]
#     plot!(figbias, [x_start, x_end], [y_start, y_end], color=:black, linewidth=2, label="")
# end

# display(figbias)