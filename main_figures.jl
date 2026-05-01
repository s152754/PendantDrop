using Plots
using LaTeXStrings
using Plots.Measures
using JLD2
using OrderedCollections
using Distributions
using Roots
using NearestNeighbors

include("src/visualization/theme.jl")
include("src/visualization/plotting.jl")
include("src/visualization/schematic.jl")
include("src/io.jl")
include("src/post_dropdata.jl")
include("src/BI_functions.jl")
include("src/priors.jl")


## user settings
CASES = ["droplet1", "droplet2", "droplet3", "droplet4", "droplet5"]
CASE = CASES[5]
prior = prior5

# placeholder for additional keyword arguments
kwargs = NamedTuple()

# paths
dir_output  = joinpath(@__DIR__, "output")
dir_figures = joinpath(@__DIR__, "figures")
dir_data    = joinpath(@__DIR__, "data", CASE)
dir_px      = joinpath(@__DIR__, "data", CASE, "allbounds")


## PD schematic
fig1 = plot_PDschem_fn(dir_output)
display(fig1)
savefig(fig1, joinpath(dir_figures, "PD_schematic.pdf"))

## Ω-study
fig2 = plot_Ω_fn(dir_output)
display(fig2)
savefig(fig2, joinpath(dir_figures, "schematic_RK4stop123.pdf"))

## distortion effect
fig3 = plot_distortion_fn(dir_data)
display(fig3)
savefig(fig3, joinpath(dir_figures, "distvsundist_$(CASE).pdf"))

## lhood data
fig4 = plot_lhood_fn(dir_output, dir_data)
display(fig4)
savefig(fig4, joinpath(dir_figures, "llhooddist.pdf"))

## bounds
fig5 = plot_dropbounds_fn(dir_px)
display(fig5)
savefig(fig5, joinpath(dir_figures, "CombinedDropletProfiles.png"))
#----------------------------------
## R-apex study
fig6 = plot_Rapex_fn(dir_data)
display(fig6)
savefig(fig6, joinpath(dir_figures, "rapex_circles.pdf"))

## std data increasing
fig7 = plot_stdincrease_fn(;legend=false)
display(fig7)
savefig(fig7, joinpath(dir_figures, "std_decrease_Ndata.pdf"))

## histogram plot
fig8 = plot_hist_fn(dir_data, "chain_out_final.jld2", prior)
display(fig8)
savefig(fig8, joinpath(dir_figures, "histprior.pdf"))

## single ppd plot
fig9 = plot(xlabel=L"r \ [\mathrm{mm}]", ylabel=L"z \ [\mathrm{mm}]")
plot_ppd_fn!(fig9, dir_output, "QoI_"*CASE*".jld2", dir_data, "chain_out_final.jld2"; aspect_ratio=:equal)
display(fig9)
savefig(fig9, joinpath(dir_figures, "ppd_"*CASE*".pdf"))

## trace plot
fig10 = plot_trace_fn(dir_data, "chain_out_final.jld2")
display(fig10)
savefig(fig10, joinpath(dir_figures, "trace.pdf"))

## ppd WO cases
fig11 = plot(xlabel=L"\mathrm{Wo} \ [-]", ylabel=L"\gamma \ [\mathrm{mN}/\mathrm{m}]")
for i = 1:length(CASES)
    dir_data_i   = joinpath(@__DIR__, "data", CASES[i])
    plot_ppd_fn!(fig11, dir_output, "QoI_"*CASES[i]*".jld2", dir_data_i, "chain_out_final.jld2"; aspect_ratio=:equal, color=mypalette[i+2], alpha=0.1, label=latexstring("\\mathrm{droplet} \\ $i"))
end
display(fig11)
savefig(fig11, joinpath(dir_figures, "ppdall.pdf"))

## model evaluations
fig12 = plot(xlabel=L"\mathrm{Wo} \ [-]", ylabel=L"\gamma \ [\mathrm{mN}/\mathrm{m}]", xlims=(-0.01, 1.01))
for i = 1:length(CASES)
    dir_data_i   = joinpath(@__DIR__, "data", CASES[i])
    plot_modelevals_fn!(fig12, dir_output, "QoI_"*CASES[i]*".jld2", dir_data_i, "chain_out_final.jld2"; color=mypalette[i+2], alpha=0.1, label=latexstring("\\mathrm{droplet} \\ $i"))
end
display(fig12)
savefig(fig12, joinpath(dir_figures, "WO_1000samples.pdf"))

## posterior WO cases
fig13 = plot(xlabel=L"\mathrm{Wo} \ [-]", ylabel=L"\gamma \ [\mathrm{mN}/\mathrm{m}]", xlims=(-0.01, 1.01))
for i = 1:length(CASES)
    dir_data_i   = joinpath(@__DIR__, "data", CASES[i])
    plot_postWO_fn!(fig13, dir_output, "QoI_"*CASES[i]*".jld2", dir_data_i, "chain_out_final.jld2"; color=mypalette[i+2], label=latexstring("\\mathrm{droplet} \\ $i"))
end
display(fig13)
savefig(fig13, joinpath(dir_figures, "WO_err.pdf"))