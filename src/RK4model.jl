using Plots, JLD2, LinearAlgebra

function ode_system(s, y, dp0, sigma, rho, g)
    r, z, psi = y[1], y[2], y[3]
    r_eff = max(r, 1e-12)
    mean_curvature = (dp0 - rho * g * z) / ( 2 * sigma)
    drds = cos(psi)
    dzds = sin(psi)
    dpsids = 2 * mean_curvature - sin(psi) / r_eff
    # return @view([drds, dzds, dpsids])[:] # small vector
    return [drds, dzds, dpsids]
end

# single RK4 step
function rk4_step(s, y, h, dp0, sigma, rho, g)
    k1 = ode_system(s, y, dp0, sigma, rho, g)
    k2 = ode_system(s + 0.5*h, y .+ 0.5*h .* k1, dp0, sigma, rho, g)
    k3 = ode_system(s + 0.5*h, y .+ 0.5*h .* k2, dp0, sigma, rho, g)
    k4 = ode_system(s + h,       y .+ h .* k3, dp0, sigma, rho, g)
    return y .+ (h/6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
end

function PD_fn(
    params::Vector{Float64}, 
    Rneedle::Float64; n_points::Int = 5000,
    rho::Float64 = 1e-3, g::Float64 = 9.81e3,
    s_max::Float64 = 20.0,
    s0::Float64 = 1e-9,
    r_tol::Float64 = 1e-12,
    h_tol::Float64 = 1e-12,
    max_bisect_iter::Int = 60,
    hits_max::Int = 2,
    Δs_min::Float64 = 1e-3,
    r_break::Float64 = 2.
    )

    sigma = exp(params[1])
    R0 = exp(params[2])
    kappa0 = 2 / R0

    dp0 =  kappa0 * sigma
    r_init = s0
    z_init = 0.0
    psi_init = kappa0 * s0
    y0 = [r_init, z_init, psi_init]

    h = (s_max - s0) / max(1, n_points - 1)

    s_vals = Float64[s0]
    y_list = Vector{Vector{Float64}}()
    push!(y_list, copy(y0))

    hits = 0
    hit_s = Float64[]

    for _ = 1:(n_points - 1)
        s_i = s_vals[end]
        y_i = y_list[end]
        r_i = y_i[1]

        # full step
        y_next = rk4_step(s_i, y_i, h, dp0, sigma, rho, g)
        r_next = y_next[1]

        # --- STOP if full-step radius exceeds limit ---
        if r_next > r_break
            @info "Stopping: r exceeded limit" r=r_next s=s_i
            if isempty(y_list)
                return Float64[], Float64[], Float64[], hits
            else
                Ypart = hcat(y_list...)'
                rpart = Ypart[:,1]
                zpart = Ypart[:,2] .- Ypart[end,2]
                spart = s_vals
                return rpart, zpart, spart, hits
            end
        end

        # crossing test
        f_i    = r_i - Rneedle
        f_next = r_next - Rneedle

        # no crossing → just accept step
        if f_i * f_next > 0
            push!(s_vals, s_i + h)
            push!(y_list, y_next)
            continue
        end

        # crossing detected -> bisection on step length
        h_left  = 0.0
        h_right = h
        h_mid   = 0.0
        y_mid   = copy(y_i)

        for _ = 1:max_bisect_iter
            h_mid = 0.5*(h_left + h_right)
            y_mid = rk4_step(s_i, y_i, h_mid, dp0, sigma, rho, g)
            f_mid = y_mid[1] - Rneedle

            # STOP if bisection radius exceeds limit
            if y_mid[1] > r_break
                @info "Stopping during bisection: r exceeded limit" r=y_mid[1] s=(s_i + h_mid)
                if isempty(y_list)
                    return Float64[], Float64[], Float64[], hits
                else
                    Ypart = hcat(y_list...)'
                    rpart = Ypart[:,1]
                    zpart = Ypart[:,2] .- Ypart[end,2]
                    spart = s_vals
                    return rpart, zpart, spart, hits
                end
            end

            if abs(f_mid) <= r_tol || (h_right - h_left) <= h_tol
                break
            end

            if f_i * f_mid < 0
                h_right = h_mid
            else
                h_left = h_mid
                f_i = f_mid
            end
        end

        # accept the crossing point
        push!(s_vals, s_i + h_mid)
        push!(y_list, copy(y_mid))
        s_hit = s_i + h_mid

        if length(hit_s) == 0
            push!(hit_s, s_hit)
            hits += 1
        elseif abs(s_hit - hit_s[end]) >= Δs_min
            push!(hit_s, s_hit)
            hits += 1
        end

        if hits == hits_max
            break
        end
    end

    # convert to arrays
    if isempty(y_list)
        return Float64[], Float64[], Float64[], hits
    end
    Y = hcat(y_list...)'   # (N,3)
    r = Y[:, 1]
    z = Y[:, 2]
    s = s_vals

    z_shifted = z .- z[end]

    return r, z_shifted, s, hits
end

# function trim_trailing_repeats(r::Vector{T}, z::Vector{S}; keepmissing::Bool = false) where {T,S}
#     length(r) == length(z) || throw(ArgumentError("r and z must have the same length"))

#     n = length(r)
#     if n == 0
#         return r, z, nothing
#     end

#     # Find the last non-missing value in r
#     lastidx = nothing
#     for i in n:-1:1
#         if keepmissing ? true : r[i] !== missing
#             lastidx = i
#             if r[i] === missing && !keepmissing
#                 continue
#             end
#             break
#         end
#     end

#     if lastidx === nothing || r[lastidx] === missing
#         return r, z, nothing
#     end

#     lastval = r[lastidx]

#     # Walk upwards from the end while r[i] == lastval
#     cut = n
#     while cut >= 1
#         if r[cut] === lastval
#             cut -= 1
#         else
#             break
#         end
#     end

#     # Keep one repeated value (cut+1)
#     r2 = r[1:cut+1]
#     z2s = z[1:cut+1]

#     # Shift z so last value becomes 0
#     z2 = z2s .- z2s[end]

#     return r2, z2, lastval
# end


# # real data
# Lp = load("data/droplet5/dropletcoords_rhs.jld2")
# rdf = Lp["rr_sorted"]
# zdf = Lp["Zr_sorted"]

# rd, zd, lastval = trim_trailing_repeats(rdf, zdf; keepmissing=false)

θs = log.([72, 1.4])
Rneedle = 0.45 #rd[end]

r, z, s, hits = PD_fn(θs, Rneedle, r_tol=1e-12, h_tol=1e-12, hits_max=2)

plot(r, z, aspect_ratio=:equal)
