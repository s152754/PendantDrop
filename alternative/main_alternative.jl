using Plots

include("src/PD_model.jl")

r, z, volume, area, kappas, kappap =PD_fn([4., 16.], rneedle=1.4, deltarho=1.2, grav=1.1) # zoals in struct

# Plot results
shape_plt = plot_shape(r, z);
curv_plt = plot_curvature(z, kappas, kappap);
combined_plt = plot(shape_plt, curv_plt, layout=(1,2); show=true)