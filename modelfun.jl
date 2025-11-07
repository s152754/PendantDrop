using Plots
using LinearAlgebra
using Interpolations
using Distributions
using Random
using JLD2

Random.seed!(1234)  # for reproducibility


include("src/PD_functions.jl")
include("src/PD_structs.jl")
include("src/utils.jl")
include("plotting.jl")


# Assign structs containing the parameters you would want to solve for
# params_phys = ParamsPhys(sigma=0.07, grav=9.81, rneedle=0.001, volume0=16e-9, deltarho=1260) # example Aricia 
sigma = 72#1*72.0          # Surface tension (mN/m)
grav = 9.81e3#4*9.81e3        # Gravity (mm/s^2)
rneedle = 0.2#2*0.3      # Needle radius (mm)
volume0 = 4.0    # Prescribed volume (mm^3)
deltarho = 1.26e-3#1.26e-3 # Density difference (rho_top - rho_bottom) (g/mm^3)
N = 100              # Number of gridpoints (numeric parameter)

params_phys = ParamsPhys(sigma=sigma, grav=grav, rneedle=rneedle, volume0=volume0, deltarho=deltarho)
# params_phys = ParamsPhys() 
params_num = ParamsNum(N=N, nMaxIter=100, epsilon=1e-6)
wo = (deltarho*grav*volume0) / (2 * π * sigma * rneedle)
println("wo = $(wo)")

# Solve for the droplet shape (Young-Laplace)
vars_sol, vars_num = gen_single_drop(params_phys, params_num; verbose=true)

# Post processing and plotting
volume, area = calculate_volume_area(vars_sol, vars_num; verbose=true);
kappas, kappap = find_curvature(vars_sol, vars_num);

# Plot results
shape_plt = plot_shape(vars_sol.r, vars_sol.z);
curv_plt = plot_curvature(vars_sol.z, kappas, kappap);
combined_plt = plot(shape_plt, curv_plt, layout=(1,2); show=true)

####################
### virtual data ###
####################
r = vars_sol.r
z = vars_sol.z
rd = rand(Normal(0,0.02*rneedle), length(r)) .+ r
zd = rand(Normal(0,0.02*rneedle), length(z)) .+ z

###################
### actual data ###
###################
L = load("data/droplet_coords_rhs_sorted_big.jld2")
rds = L["rr_sorted"]
zds = L["Zr_sorted"]

## estimate the volume ##
itp = LinearInterpolation(zds, rds, extrapolation_bc=Flat())

# Define integration range
zmin = minimum(zds)
zmax = maximum(zds)

# Discretize for numerical integration
z_vals = range(zmin, zmax, length=1000)
r_vals = itp.(z_vals)

# Compute volume using trapezoidal rule
volume = π * sum(diff(z_vals) .* ((r_vals[1:end-1].^2 + r_vals[2:end].^2) ./ 2))

println("Droplet volume ≈ $(volume) mm³")



fac = 1.0 #rneedle / rds[end]
rd = rds .* fac
zd = zds .* fac

facfig = 1.2
figcomp = plot(xlabel="r", ylabel="z")
ymax = maximum([maximum(zd)*facfig, maximum(z)*facfig])
ymin = minimum([minimum(zd)*facfig, minimum(z)*facfig])
xmax = maximum([maximum(rd)*facfig, maximum(r)*facfig])
xmin = minimum([minimum(rd)*facfig, minimum(r)*facfig])
plot!(figcomp, rd, zd; markershape=:circle, grid=true, xlim=[xmin, xmax], ylim=[ymin, ymax], aspect_ratio=:equal, label="data") # xlim=[-0.01, maximum(r)*1.1]
plot!(figcomp, r, z; markershape=:circle, linewidth=2, label="model")
plot!(figcomp, legend=:topright)
# plot!(figcomp, rd, zd; markershape=:circle, grid=true, aspect_ratio=:equal) # xlim=[-0.01, maximum(r)*1.1]
display(figcomp)