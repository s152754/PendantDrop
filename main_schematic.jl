using JLD2

include("src/RK4model.jl")

# PD schematic
Rneedle = 0.45
Rapex = 1.4
γ = 72.0

θs = [γ, Rapex]
r, z, s, hits = PD_fn(log.(θs), Rneedle, n_points=10_000, r_tol=1e-12, h_tol=1e-12, hits_max=2)
jldsave("output/schematic.jld2"; r, z, s, hits, Rapex, γ, Rneedle)


# Ω - study
Rneedle = 0.45
Ωs      = 3
Rapex   = 1.2
γ       = 72.0
θs      = [γ, Rapex]

for i = 1:Ωs
    r, z, s, hits = PD_fn(log.(θs), Rneedle, n_points=5000, r_tol=1e-12, h_tol=1e-12, hits_max=i)
    jldsave("output/schematic_Ω$(hits).jld2"; r, z, s, hits, Rapex, γ, Rneedle)
end