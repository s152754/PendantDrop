using JLD2, Plots, LaTeXStrings, OrderedCollections, Distributions
using NearestNeighbors
using Plots.Measures

# pgfplotsx()  # switch backend

# using PGFPlotsX

paththesis = "C:\\Users\\s152754\\VScode\\PhDproject\\adjustthesis\\chapter5\\"
mypalette = ["#882255", "#661100", "#6699CC", "#888888", "#332288", "#AA4499", "#44AA99", "#999933", "#88CCEE", "#CC6677", "#DDCC77", "#117733"]

include("src/post_dropdata.jl")
include("src/io.jl")
include("src/BI_functions.jl")
include("src/visualization/plotting.jl")

case = 5
Lp = load("data/droplet$(case)/dropletcoords_out.jld2")
rdf = Lp["rs_sorted"]
zdf = Lp["Zs_sorted"]

Lpdist = load("data/droplet$(case)/dropletcoords_out_dist.jld2")
rdf_dist = Lpdist["rs_sorted"]
zdf_dist = Lpdist["Zs_sorted"]

rd_undist, zd_undist    = reorder_droplet(rdf, zdf)
rd_dist, zd_dist        = reorder_droplet(rdf_dist, zdf_dist)
test = distortion_fn(hcat([rd_undist, zd_undist]...), hcat([rd_dist, zd_dist]...))
display(test)
# rd_dist2, zd_dist2    = trim_top_duplicates(rd_dist, zd_dist)


figori = plot(xlabel=L"r \ [\mathrm{mm}]", ylabel=L"z \ [\mathrm{mm}]", aspect_ratio=:equal)
plot!(legend=:topright, legend_font_halign=:left)
# plot!(figori, size=(600, 410), labelfontsize=23, tickfontsize=17, legendfontsize=19, bottom_margin=1mm)
# plot!(figori, size=(485, 400), labelfontsize=25, tickfontsize=19, legendfontsize=21, xlims=(-1.6,3.4), ylims=(-4,0), aspect_ratio=1)
plot!(figori, rd_dist, zd_dist, lw=2, label=L"\mathrm{Distorted}", color=mypalette[10])
plot!(figori, rd_undist, zd_undist, lw=2, label=L"\mathrm{Undistorted}", color=mypalette[7])
display(figori)

# ms = hcat([rd_undist, zd_undist]...)  
# idxs, deltas = find_nearest_neighbors(ms, rd_dist, zd_dist)

# rmse = sqrt(mean(hcat(deltas...).^2))

# pathout = "C:\\Users\\s152754\\OneDrive - TU Eindhoven\\Documents\\PhD\\Project\\Publications\\Paper 4\\figures\\"
# # savefig(figori, pathout*"distvsundist_droplet5.pdf")
# savefig(figori, paththesis*"distvsundist_droplet5_test.pdf")