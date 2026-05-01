# distortion effect
function plot_distortion_fn(data_dir; kwargs...)

    isdir(data_dir) || error("Directory not found: $data_dir")

    # --- load undistorted ---
    file_undist = joinpath(data_dir, "dropletcoords_out.jld2")
    @load file_undist rs_sorted Zs_sorted

    rdf = rs_sorted
    zdf = Zs_sorted

    # --- load distorted ---
    file_dist = joinpath(data_dir, "dropletcoords_out_dist.jld2")
    @load file_dist rs_sorted Zs_sorted

    rdf_dist = rs_sorted
    zdf_dist = Zs_sorted

    # --- reorder ---
    rd_undist, zd_undist = reorder_droplet(rdf, zdf)
    rd_dist,   zd_dist   = reorder_droplet(rdf_dist, zdf_dist)

    fig = plot(xlabel=L"r \ [\mathrm{mm}]", 
            ylabel=L"z \ [\mathrm{mm}]", 
            aspect_ratio=:equal; 
            kwargs...)

    plot!(fig, rd_undist, zd_undist, label=L"\mathrm{Undistorted}", color=mypalette[7], markersize=3)
    plot!(fig, rd_dist, zd_dist, label=L"\mathrm{Distorted}", color=mypalette[10], markersize=3)  

    return fig
end

## trace plot
function plot_trace_fn(data_dir, filename; kwargs...)

    isdir(data_dir) || error("Directory not found: $data_dir")

    file = joinpath(data_dir, filename)
    @load file flatchain chain

    flatchainlin = exp.(flatchain)
    numdims = size(flatchain)[1]
    n_cols = ceil(Int, sqrt(numdims))
    n_rows = ceil(Int, numdims / n_cols)

    labels = [L"\gamma \ [\mathrm{mN}/\mathrm{m}]", L"R_\mathrm{apex} \ [\mathrm{mm}]"]

    fig = plot(layout=(n_rows, n_cols), xlabel=L"N / K  \ [-]", legend=false; kwargs...)
    for i in eachindex(labels)
        plot!(fig, exp.(chain[i,:,:])', ylabel=labels[i], subplot=i, left_margin=10mm)
    end

    return fig
end


function plot_hist_fn(data_dir, filename, prior; kwargs...)

    isdir(data_dir) || error("Directory not found: $data_dir")

    file = joinpath(data_dir, filename)
    @load file flatchain chain

    flatchainlin = exp.(flatchain)
    numdims = size(flatchain)[1]
    n_cols = ceil(Int, sqrt(numdims))
    n_rows = ceil(Int, numdims / n_cols)

    labels = [L"\gamma \ [\mathrm{mN}/\mathrm{m}]", L"R_\mathrm{apex} \ [\mathrm{mm}]"]

    fig = plot(layout=(n_rows, n_cols), legend=false; kwargs...)

    for i in eachindex(labels)
            xtick_positions = range(minimum(flatchainlin[i,:]), stop=maximum(flatchainlin[i,:]), length=3)
            xtick_labels = [string(round(x, sigdigits=3)) for x in xtick_positions]
            histogram!(fig, flatchainlin[i,:], xlabel=labels[i], subplot=i, normalize=:pdf, color=:black, xticks = (xtick_positions, xtick_labels), label="")
            xp = range(minimum(flatchain[i,:]), stop=maximum(flatchain[i,:]), length=200)
            xplin = exp.(xp)
            J = exp.(xp)
            plot!(fig, xplin, pdf.(collect(values(prior))[i], xp)./J, subplot=i,color=:red, lw=2, label="")
    end

    return fig
end


function plot_ppd_fn!(fig, dir_output, filename_out, dir_data, filename_chain; kwargs...)

    # load
    file_out = joinpath(dir_output, filename_out)
    @load file_out QoIppd 

    file_chain = joinpath(dir_data, filename_chain)
    @load file_chain flatchain

    file_data = joinpath(dir_data, "dropletcoords_rhs_out.jld2")
    @load file_data rr_sorted Zr_sorted

    rdfull, zdfull, lastval = trim_trailing_repeats(rr_sorted, Zr_sorted; keepmissing=false)

    Nrep = 16
    ds, sdfull = determine_arclength(rdfull, zdfull)
    rd, zd = trim_data(Nrep, rdfull, zdfull, sdfull)

    Nsamples = size(QoIppd)[1]

    samplesppd = Vector{Matrix{Float64}}(undef, Nsamples)
    x_mean, y_mean, x_upper, y_upper, x_lower, y_lower = QoIppd_plotting(Nsamples, samplesppd, QoIppd, flatchain)

    band_x = vcat(x_upper, reverse(x_lower))
    band_y = vcat(y_upper, reverse(y_lower))
    plot!(fig, band_x, band_y, seriestype=:shape, color=mypalette[7], alpha=0.5; kwargs...)
    scatter!(fig, rd, zd, markersize=2, color=:black)

    return fig
end


function plot_stdincrease_fn(; kwargs...)

    Ndata = [550, 82, 55, 28, 16]
    σdatap = [0.189, 0.520, 0.644, 0.959, 1.23]

    fig = plot(xaxis=:log, yaxis=:log, xlabel=L"N_\mathrm{data} \ [-]", ylabel=L"\sigma_{\mathrm{post}, \gamma} \ [\mathrm{mN}/\mathrm{m}]"; kwargs...)
    plot!(fig, Ndata, σdatap, marker=:circle, markersize=5, lw=2, color=:black)
    plot!(fig, Ndata, 7.0 .* Ndata.^(-0.5), lw=2, color=:black, linestyle=:dash)
    annotate!(fig, 110, 0.8, text(latexstring("-\\frac{1}{2}"), 20, :left))

    return fig
end


function plot_modelevals_fn!(fig, dir_output, filename_out, dir_data, filename_chain; Nwo=1000, kwargs...)
    
    # load
    file_out = joinpath(dir_output, filename_out)
    @load file_out Wo

    file_chain = joinpath(dir_data, filename_chain)
    @load file_chain chain flatchain

    flatchainlin = exp.(flatchain)

    scatter!(fig, Wo[end-Nwo+1:end], flatchainlin[1,end-Nwo+1:end]; kwargs...)
    
    return fig
end


function plot_postWO_fn!(fig, dir_output, filename_out, dir_data, filename_chain; Nwo=1000, kwargs...)
    
    # load
    file_out = joinpath(dir_output, filename_out)
    @load file_out μWo

    file_chain = joinpath(dir_data, filename_chain)
    @load file_chain chain flatchain

    flatchainlin = exp.(flatchain)
    μγ = mean(flatchainlin[1,:])
    σγ = std(flatchainlin[1,:])

    scatter!(fig, [μWo], [μγ], yerr=[2*σγ]; kwargs...)
    
    return fig
end