using Plots, JLD2, LinearAlgebra, Distributions, OrderedCollections, NearestNeighbors
using Statistics
using Roots
using LaTeXStrings
using Plots.Measures
using DataFrames
using CSV

include("postprocessing_fns.jl")
include("src/visualization/plot_helpers.jl")
include("src/RK4model.jl")

# pgfplotsx()  # switch backend

# using PGFPlotsX


paththesis = "C:\\Users\\s152754\\VScode\\PhDproject\\adjustthesis\\chapter5\\"

function arrow0!(x, y, u, v; as=0.07, lc=:black, la=1)
    nuv = sqrt(u^2 + v^2)
    v1, v2 = [u;v] / nuv,  [-v;u] / nuv
    v4 = (3*v1 + v2)/3.1623  # sqrt(10) to get unit vector
    v5 = v4 - 2*(v4'*v2)*v2
    v4, v5 = as*nuv*v4, as*nuv*v5
    plot!([x,x+u], [y,y+v], lc=lc,la=la)
    plot!([x+u,x+u-v5[1]], [y+v,y+v-v5[2]], lc=lc, la=la)
    plot!([x+u,x+u-v4[1]], [y+v,y+v-v4[2]], lc=lc, la=la)
end

function mycircle!(plt, cx, cy, r;
                        start_angle::Real=0.0,
                        end_angle::Real=2π,
                        n::Int=200,
                        color=:black,
                        lw::Real=1.5,
                        )
    θ = range(start_angle, end_angle; length=n)
    xs = cx .+ r .* cos.(θ)
    ys = cy .+ r .* sin.(θ)
    plot!(plt, xs, ys; color=color, lw=lw, label=latexstring("R_\\mathrm{apex} = $(r) \\ \\mathrm{mm}")  )
    return plt
end

mypalette = ["#882255", "#661100", "#6699CC", "#888888", "#332288", "#AA4499", "#44AA99", "#999933", "#88CCEE", "#CC6677", "#DDCC77", "#117733"]

################################
### include relevant scripts ###
################################

include("src/RK4model.jl")
include("src/BI_functions.jl")
include("src/io.jl")
include("src/priors.jl")
include("src/calculate_area_volume.jl")
include("src/post_dropdata.jl")


pathout = "C:\\Users\\s152754\\OneDrive - TU Eindhoven\\Documents\\PhD\\Project\\Publications\\Paper 4\\figures\\"
global hits = 0
######################
### single droplet ###
######################
prior = prior5

# exp. data
Lp = load("data/droplet5/dropletcoords_rhs_out.jld2")
rdf = Lp["rr_sorted"]
zdf = Lp["Zr_sorted"]

rdfull, zdfull, lastval = trim_trailing_repeats(rdf, zdf; keepmissing=false)

Nrep        = 16
ds, sdfull  = determine_arclength(rdfull, zdfull)
rd, zd, sd  = trim_data(Nrep, rdfull, zdfull, sdfull)
Rneedle     = rd[end]

# ## figure arclength S
# figs = plot(xlabel=L"N_{\mathrm{data}} \ [-]", ylabel=L"s_{\mathrm{norm}} \ [-]", legend=:topright)
# Nreps = [16, 28, 55, 82, 550]
# for i = 1:length(Nreps)
#     Nrep = Nreps[i]
#     rd, zd, sd = trim_data(Nrep, rdfull, zdfull, sdfull)
#     s_norm = (sd[end]) ./ sdfull[end]
#     scatter!(figs, [Nrep], [s_norm]; label="N_{\\mathrm{data}} = $(Nrep)", markersize=4)
# end
# display(figs)
# Rneedle = rd[end]

## chains
L = load("data/droplet5/chain_out_final.jld2")

flatchain = L["flatchain"]
chain = L["chain"]

df = DataFrame(Matrix(exp.(flatchain)'), :auto)
CSV.write("C:\\Users\\s152754\\PycharmProjects\\nutils-squeezeflow\\cornerdata_droplet5_final.csv", df)

modelfn         = θs -> PD_fn(θs, Rneedle, n_points=5000, r_tol=1e-12, h_tol=1e-12, hits_max=2)# returns [r,z]

# plotting

labels = [L"\gamma \ [\mathrm{mN}/\mathrm{m}]", L"R_\mathrm{apex} \ [\mathrm{mm}]"]

flatchainlin = exp.(flatchain)
numdims = size(flatchain)[1]
n_cols = ceil(Int, sqrt(numdims))
n_rows = ceil(Int, numdims / n_cols)

## trace plot
pt = plot(layout=(n_rows, n_cols), xlabel=L"N / K  \ [-]", legend=false)
# plot!(pt, guidefont=18, tickfont=14, legendfont=14, size=(800, 400), left_margin=5mm, right_margin=5mm, bottom_margin=6mm)
plot!(pt, labelfontsize=18, tickfontsize=14, legendfontsize=15)

for i in eachindex(labels)
    plot!(pt, exp.(chain[i,:,:])', ylabel=labels[i], subplot=i, left_margin=10mm)
end
display(pt)
# savefig(pt, pathout*"trace.pdf")
# savefig(pt, paththesis*"trace_test.pdf")

# # histogram vs prior
# p = plot(layout=(n_rows, n_cols), legend=false, size=(600, 400), guidefont=18, tickfont=14, legendfont=14, right_margin=2mm, bottom_margin=2mm)

# for i in eachindex(labels)
#         xtick_positions = range(minimum(flatchainlin[i,:]), stop=maximum(flatchainlin[i,:]), length=3)
#         xtick_labels = [string(round(x, sigdigits=3)) for x in xtick_positions]
#         histogram!(p, flatchainlin[i,:], xlabel=labels[i], subplot=i, normalize=:pdf, color=:black, xticks = (xtick_positions, xtick_labels), label="")
#         xp = range(minimum(flatchain[i,:]), stop=maximum(flatchain[i,:]), length=200)
#         xplin = exp.(xp)
#         J = exp.(xp)
#         plot!(p, xplin, pdf.(collect(values(prior))[i], xp)./J, subplot=i,color=:red, lw=2, label="")
# end
# display(p)   
# # savefig(p, pathout*"histprior.pdf")


# # ppd plot
ptl         = 0.9/139 # 128
σn          = 3*ptl
Σ           =  σn
Nsamples    = 2000
QoI, QoIppd      = compute_ppd(Nsamples, modelfn, flatchain, rd, zd, Σ)

samplesppd = Vector{Matrix{Float64}}(undef, Nsamples)
x_mean, y_mean, x_upper, y_upper, x_lower, y_lower = QoIppd_plotting(Nsamples, samplesppd, QoIppd, flatchain)

# figQoI = plot(xlabel=L"r \ [\mathrm{mm}]", ylabel=L"z \ [\mathrm{mm}]", aspect_ratio=:equal)
# plot!(figQoI, rd, zd)
# for i = 1:Nsamples
#     plot!(figQoI, QoI[i][:,1], QoI[i][:,2], color=:red, lw=3, label="", linewidth=2, legend=:topright)
# end
# display(figQoI)
##
figppdplot = plot(
    x_mean,
    y_mean,
    color=mypalette[7],
    legend=false,
    linewidth=1,
    xlabel=L"r \, [\mathrm{mm}]",
    ylabel=L"z \, [\mathrm{mm}]",
    labelfontsize=20,
    tickfontsize=14,
    size=(1.14*200, 400),
    xlims=(0,2),
    ylims=(-4.1,0.1)
)
# plot!(figppdplot, size=(300, 400), aspect_ratio=:equal, guidefont=18, tickfont=14, legendfont=14, left_margin=-10mm, bottom_margin=-3mm)

# plot!(figppdplot, size=(1.22*300, 600), xlims=(-0.1,2), ylims=(-4.1,0)) # 430, 800
# plot!(figppdplot, xlims=(-0.05,2), ylims=(-3.8, 0.05))
band_x = vcat(x_upper, reverse(x_lower))
band_y = vcat(y_upper, reverse(y_lower))
plot!(figppdplot, band_x, band_y, seriestype=:shape, color=mypalette[7], alpha=0.5)
scatter!(figppdplot, rd, zd, markersize=2, color=:black)
display(figppdplot)
#  savefig(figppdplot, pathout*"pdd_drop5.pdf")
# savefig(figppdplot, paththesis*"pdd_drop5_test.pdf")

##
####################
### all droplets ###
####################

# #####################
# ### rapex circles ###
#####################
Lptot = load("data/droplet5/dropletcoords_out.jld2")
rptot = Lptot["rs_sorted"]
zptot = Lptot["Zs_sorted"]

rptots, zptots = reorder_droplet(rptot, zptot)
rptotf, zptotf = trim_top_duplicates(rptots, zptots)

Rapexs = 2*[0.4, 0.7, 1.0]
figrapex = plot(xlabel=L"r \, [\mathrm{mm}]", ylabel=L"z \, [\mathrm{mm}]", legend=:topright, legend_font_halign=:left)
plot!(figrapex, size=(550, 415), guidefont=20, tickfont=14, legendfont=15)
plot!(figrapex, rptotf, zptotf, label=L"\mathrm{Experimental} \ \mathrm{data}", color=:black, lw=2)
plot!(figrapex, xlims=(-1.5,4))
for i = 1:3
    Rapex = Rapexs[i]
    mycircle!(figrapex, 0.0, minimum(zptotf)+Rapex, Rapex, start_angle=1.25π, end_angle=1.75π; color=mypalette[i+4], lw=2)#latexstring("R_\\mathrm{apex} = $(Rapex) \ \\mathrm{mm}"))
end
display(figrapex)
# savefig(figrapex, pathout*"rapex_circles.pdf")
# savefig(figrapex, paththesis*"rapex_circles_test.pdf")

##
#########################################
### std decrease thinning data points ###
#########################################
Ndata = [550, 82, 55, 28, 16]
σdatap = [0.189, 0.520, 0.644, 0.959, 1.23]

figstddata = plot(xaxis=:log, yaxis=:log, xlabel=L"N_\mathrm{data} \ [-]", ylabel=L"\sigma_{\mathrm{post}, \gamma} \ [\mathrm{mN}/\mathrm{m}]")
# plot!(figstddata, legend=false, size=(600, 400), guidefont=18, tickfont=14, legendfont=14, bottom_margin=3mm, left_margin=2mm)
plot!(figstddata, legend=false, size=(600, 400), labelfontsize=20, tickfontsize=16, legendfontsize=17)
plot!(figstddata, Ndata, σdatap, marker=:circle, markersize=5, lw=2, color=:black)
plot!(figstddata, Ndata, 7.0 .* Ndata.^(-0.5), lw=2, color=:black, linestyle=:dash)
annotate!(figstddata, 110, 0.8, text(latexstring("-\\frac{1}{2}"), 20, :left))
display(figstddata)
# savefig(figstddata, pathout*"std_decrease_Ndata.pdf")
# savefig(figstddata, paththesis*"std_decrease_Ndata_test.pdf")
##

##################################################
### create filled droplet from rhs coordinates ###
##################################################
Nwo = 1000
WOS = []
FLATS = []
figwo = plot(xlabel=L"\mathrm{Wo} \ [-]", ylabel=L"\gamma \ [\mathrm{mN}/\mathrm{m}]")
# plot!(figwo, size=(600, 400), guidefont=18, tickfont=14, legendfont=14, right_margin=5mm)
plot!(figwo, legend=:topright, labelfontsize=20, tickfontsize=16, legendfontsize=17, legend_font_halign=:left)
plot!(figwo, xlims=(-0.01, 1.01))
for j = 1:5

    # chains
    L = load("data/droplet$(j)/chain_out_final.jld2")

    flatchain = L["flatchain"]
    flatchainlin = exp.(flatchain)
    chain = L["chain"]
    Wos = []

    println("mean = $(mean(flatchainlin[1,end-Nwo+1:end])), std = $(std(flatchainlin[1,end-Nwo+1:end]))")

    for i = 1:Nwo
        rm, zm = modelfn(vec(flatchain[:,i]))
        coords = [(rm[j], zm[j]) for j in 1:size(rm)[1]]
        volume = axisymmetric_volume_right(coords)
        wo = wo_fn(volume, flatchainlin[1,end-i])
        # if i ==1 || i == 10
        #     println("volume = $(volume)", "γ = $(flatchainlin[1,end-i])", "wo = $(wo)")
        # end
    
        push!(Wos, wo)
    end

    scatter!(figwo, Wos, flatchainlin[1,end-Nwo+1:end], alpha=0.1, color=mypalette[j+2], label="")
    # scatter!([NaN], [NaN], label=latexstring("\\mathrm{droplet} \\ $(j)"), color=mypalette[j+2])
    push!(WOS, Wos)
    push!(FLATS, flatchainlin[1,end-Nwo+1:end])
end

scatter!(figwo, [2], [FLATS[1][1]], alpha=1, color=mypalette[3], label=L"\mathrm{Droplet} \ 5")
scatter!(figwo, [2], [FLATS[2][1]], alpha=1, color=mypalette[4], label=L"\mathrm{Droplet} \ 4")
scatter!(figwo, [2], [FLATS[3][1]], alpha=1, color=mypalette[5], label=L"\mathrm{Droplet} \ 3")
scatter!(figwo, [2], [FLATS[4][1]], alpha=1, color=mypalette[6], label=L"\mathrm{Droplet} \ 2")
scatter!(figwo, [2], [FLATS[5][1]], alpha=1, color=mypalette[7], label=L"\mathrm{Droplet} \ 1")

display(figwo)
# savefig(figwo, pathout*"WO_1000samples.pdf")
# savefig(figwo, paththesis*"WO_1000samples_test.pdf")

## WO via exp

figwoerr = plot(xlabel=L"\mathrm{Wo} \ [-]", ylabel=L"\gamma \ [\mathrm{Pa} \cdot \mathrm{s}]")
# plot!(figwoerr, size=(600, 400), guidefont=18, tickfont=14, legendfont=14, right_margin=5mm, legend=:bottomright)
plot!(figwoerr, labelfontsize=20, tickfontsize=16, legendfontsize=17, legend=:bottomright, legend_font_halign=:left)
plot!(figwoerr, xlims=(-0.01, 1.01))
for j = 1:5

    # exp. data
    Lp = load("data/droplet$(j)/dropletcoords_rhs_out.jld2")
    rdf = Lp["rr_sorted"]
    zdf = Lp["Zr_sorted"]

    rdfull, zdfull, lastval = trim_trailing_repeats(rdf, zdf; keepmissing=false)

    Nrep = 16
    ds, sdfull = determine_arclength(rdfull, zdfull)
    rd, zd = trim_data(Nrep, rdfull, zdfull, sdfull)

    Rneedle = rd[end]

    

    # chains
    L = load("data/droplet$(j)/chain_out_final.jld2")

    flatchain = L["flatchain"]
    flatchainlin = exp.(flatchain)

    μγ = mean(flatchainlin[1,:])
    σγ = std(flatchainlin[1,:])


    coords = [(rd[j], zd[j]) for j in 1:size(rd)[1]]
    volume = axisymmetric_volume_right(coords)

    wo = wo_fn(volume, μγ)
    
    # plot!(figwoerr, [wo], [μγ], yerr=2*σγ)
    
    scatter!(figwoerr, [wo], [μγ], yerr=2*σγ; marker=:circle, markersize=5, color=mypalette[j+2], label=latexstring("\\mathrm{droplet} \\ $(j)"))


end
display(figwoerr)
# savefig(figwoerr, pathout*"WO_err.pdf")
# savefig(figwoerr, paththesis*"WO_err_test.pdf")


##############################
### compute QoI and QoIppd ###
##############################
# global QoI = []
# global QoIppd = []




figpddplot = plot(xlabel="r [mm]", ylabel="z [mm]", legend=:left)

CVs = []
Volumes = []
WOs = []
μγs = []
σγs = []
priors = [prior1, prior2, prior3, prior4, prior5]
cases = [1, 2, 3, 4, 5]
figproj= plot(xlabel="r [mm]", ylabel="z [mm]", aspect_ratio=:equal, legend=:topleft)
plot!(figproj, size=(600, 400))
# for i = 1:5
i = 5
QoI = []
QoIppd = []

# exp. data
Lp = load("data/droplet$(i)/dropletcoords_rhs_out.jld2")
rdf = Lp["rr_sorted"]
zdf = Lp["Zr_sorted"]

rdfull, zdfull, lastval = trim_trailing_repeats(rdf, zdf; keepmissing=false)

Nrep = 16
ds, sdfull = determine_arclength(rdfull, zdfull)
rd, zd = trim_data(Nrep, rdfull, zdfull, sdfull)

R_cap = rd[end]

modelfn     = θs -> PD_fn(θs, R_cap, n_points=10000, r_tol=1e-12, h_tol=1e-12, hits_max=2)# returns [r,z]

# chains
L = load("data/droplet$(i)/chain_out_final.jld2")

flatchain = L["flatchain"]
chain = L["chain"]

##
##################
### llhooddist ###
##################

r_model, zm, sm, hitsm        = modelfn(log.([72.0, 1.4]))
# figllhooddist = llhooddist_fig(r_model, zm, rd, zd)
(figllhooddist,zm_zoom, index_zm, zd_zoom) = llhooddist_figgg(r_model, zm, rd, zd)
display(figllhooddist)
# savefig(final_plot, pathout*"llhooddist.pdf")
# savefig(figllhooddist, paththesis*"llhooddist_test.pdf")

##
###############
### all ppd ###
###############
figppdall = plot(xlabel=L"r \ [\mathrm{mm}]", ylabel=L"z \ [\mathrm{mm}]")
plot!(figppdall, size=(450, 400), guidefont=18, tickfont=14, legendfont=14)
plot!(figppdall, xlims=(-0.05,4.5), ylims=(-3.8, 0.05))

for i = 1:5
    j = 5:-1:1
    # exp. data
    Lp = load("data/droplet$(i)/dropletcoords_rhs_out.jld2")
    rdf = Lp["rr_sorted"]
    zdf = Lp["Zr_sorted"]

    rdfull, zdfull, lastval = trim_trailing_repeats(rdf, zdf; keepmissing=false)

    Nrep = 16
    ds, sdfull = determine_arclength(rdfull, zdfull)
    rd, zd = trim_data(Nrep, rdfull, zdfull, sdfull)

    # println("rd = $(size(rd)), zd = $(size(zd))")

    # chain
    L = load("data/droplet$(i)/chain_out_final.jld2")

    flatchain = L["flatchain"]
    chain = L["chain"]

    flatchainlin = exp.(flatchain)

    # QoI ppd
    QoI, QoIppd      = compute_ppd(Nsamples, modelfn, flatchain, rd, zd, Σ)

    samplesppd = Vector{Matrix{Float64}}(undef, Nsamples)
    x_mean, y_mean, x_upper, y_upper, x_lower, y_lower = QoIppd_plotting(Nsamples, samplesppd, QoIppd, flatchain)
    band_x = vcat(x_upper, reverse(x_lower))
    band_y = vcat(y_upper, reverse(y_lower))

    
    scatter!(figppdall, rd, zd, markershape=:circle, label="", markersize=1, color=:black)
    plot!(figppdall, x_mean, y_mean, color=mypalette[i+2], label="", linewidth=1, legend=:topright, legend_font_halign=:left)
    plot!(figppdall, band_x, band_y, seriestype=:shape, color=mypalette[i+2], alpha=0.5, label=latexstring("\$\\mathrm{Droplet} \\ $(j[i])\$"))

end
scatter!(figppdall, [-1], [-1], label=L"\mathrm{Experimental} \ \mathrm{data}", color=:black)
display(figppdall)
# savefig(figppdall, pathout*"ppdall.pdf")
# savefig(figppdall, paththesis*"ppdall_test.pdf")






  
#     ##########################
#     ### posterior vs prior ###
#     ##########################
#     prior = priors[i]

#     flatchainlin = exp.(flatchain)

#     QoI, QoIppd      = compute_ppd(Nsamples, modelfn, flatchain, rd, zd, Σ)

#     samplesppd = Vector{Matrix{Float64}}(undef, 100)
#     x_mean, y_mean, x_upper, y_upper, x_lower, y_lower = QoIppd_plotting(samplesppd, QoIppd, flatchain)

#     # Plot
#     plot!(figpddplot, x_mean, y_mean, color=:red, lw=3, label="", linewidth=2, legend=:topright)
#     band_x = vcat(x_upper, reverse(x_lower))
#     band_y = vcat(y_upper, reverse(y_lower))
#     plot!(figpddplot, aspect_ratio=:equal, xlims=(-0.01, 1.5), ylims=(-1, 0.01))
#     plot!(figpddplot, band_x, band_y, seriestype=:shape, color=:red, alpha=0.3, label="")
#     scatter!(figpddplot, rd, zd, label="", markersize=2, color=:black)
    
#     println("case $(i): mean = $(mean(flatchainlin, dims=2)), std = $(std(flatchainlin, dims=2)), CV = $(std(flatchainlin, dims=2) ./ mean(flatchainlin, dims=2) * 100)")

#     # CV based graph
#     μγ = mean(flatchainlin[1,:])
#     σγ = std(flatchainlin[1,:])
#     CVγ = σγ / μγ * 100

#     coords = [(rd[j], zd[j]) for j in 1:length(rd)]
#     volume = axisymmetric_volume_right(coords)
#     wo = wo_fn(volume, mean(flatchainlin[1,:]))
#     println("wo case $(i) = $(wo)")
#     push!(Volumes, volume)
#     push!(CVs, CVγ)
#     push!(WOs, wo)
#     push!(σγs, σγ)
#     push!(μγs, μγ)
# end

# scatter!(figpddplot, [NaN], [NaN], label="exp. data", color=:black)
# plot!(figpddplot, [NaN], [NaN], label="ppd", color=:red)
# display(figpddplot)

# figcv = plot(xlabel="V [mm3]", ylabel="CV_γ [-]")
# ratem = -1.7
# yrate = 10*Volumes.^(ratem)
# scatter!(figcv, Volumes, CVs, xaxis=:log, yaxis=:log, label="")
# plot!(figcv, Volumes, yrate, label="rate = $(ratem)")
# display(figcv)

# figwo = plot(xlabel="WO [-]", ylabel="γ [Pa*s]")
# for i = 1:5
#     scatter!(figwo, [WOs[i]], [μγs[i]], yerr=[2*σγs[i]], label="case $(i)", markersize=3)
# end

# display(figwo)
# savefig(figwo, pathout*"WO.pdf")



# # coords = [(rd[i], zd[i]) for i in 1:length(rd)]
# # include("src/calculate_area_volume.jl")

# # println("Axisymmetric Volume: ", axisymmetric_volume_right(coords))












