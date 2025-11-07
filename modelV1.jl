using Plots
using LinearAlgebra
using Interpolations


include("src/PD_functions.jl")
include("src/PD_structs.jl")
include("src/utils.jl")
include("plotting.jl")


# Assign structs containing the parameters you would want to solve for
# params_phys = ParamsPhys(sigma=72, grav=9.81e3, rneedle=1., volume0=32, deltarho=1.26e-6) # example Aricia 
params_phys = ParamsPhys() 
params_num = ParamsNum()

# Solve for the droplet shape (Young-Laplace)
vars_sol, vars_num = gen_single_drop(params_phys, params_num; verbose=true)

# Post processing and plotting
volume, area = calculate_volume_area(vars_sol, vars_num; verbose=true);
kappas, kappap = find_curvature(vars_sol, vars_num);

# Plot results
shape_plt = plot_shape(vars_sol.r, vars_sol.z);
curv_plt = plot_curvature(vars_sol.z, kappas, kappap);
combined_plt = plot(shape_plt, curv_plt, layout=(1,2); show=true)