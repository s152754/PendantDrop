include(joinpath(@__DIR__, "plot_helpers.jl"))
include(joinpath(@__DIR__, "../RK4model.jl"))


function plot_Ω_fn(data_dir; kwargs...)
    
    files = filter(f -> occursin("Ω", f), readdir(data_dir; join=true))
    sort!(files)

    fig = plot(xlabel=L"r \ [\mathrm{mm}]", 
                ylabel=L"z \ [\mathrm{mm}]",
                ylims=(-4.0,0), 
                xlims=(0,4),
                aspect_ratio=1)

    # Load Rneedle once from the first file
    @load files[1] Rneedle
    vline!(fig, [Rneedle], ls=:dash, color=:black, label=latexstring("R_\\mathrm{n}  = $(Rneedle) \\ \\mathrm{mm}"))

    for file in files
        @load file r z s hits
        plot!(fig, r, z, lw=3, color=mypalette[hits+4], label=latexstring("\\Omega = $(hits)"))
    end

    return fig
end


function plot_PDschem_fn(datadir; kwargs...)

    # load file
    file = joinpath(datadir, "schematic.jld2")
    @load file r z s hits Rapex Rneedle


    rtot = vcat(-reverse(r), r)
    ztot = vcat(reverse(z), z) .- z[1]

    si          = 70
    sf          = 350
    sline       = hcat(r[si:sf] .-0.05, z[si:sf].- z[1] .+ 0.05)
    idx         = max(2, min(1050, length(r) - 1))  # Clamp index to valid range
    tanx, tany  = tangent_to_y0(hcat(r, z .-z[1]), idx)


    fig = plot(aspect_ratio=:equal,
        xticks = false,
        yticks = false,
        grid = false,
        framestyle=:none,
        legend=false)


    # gray axis
    add_axes_with_arrows!(fig; xlims=(0.0, maximum(rtot)*1.1), ylims=(minimum(ztot), maximum(ztot)), arrow_len=0.35, color_axes=:gray, color_arrows=:black)                  

    # black axis
    add_axes_with_arrows_short!(fig; xlims=(0.0, maximum(rtot)*1.1), ylims=(minimum(ztot), maximum(ztot)), arrow_len=0.35, color=:black)   

    # s axis
    plot!(fig, sline[:,1], sline[:,2], color=:black, arrow=:closed)
    annotate!(fig, maximum(sline[:,1])+0.05, maximum(sline[:,2])+0.1, text(L"s", :black))

    # model
    plot!(fig, rtot, ztot, color=mypalette[3], lw=3)

    # needle
    dashed_ellipse!(fig, 0.0, ztot[1]; r=rtot[end], scale_y=0.2, angle=0.0, n=200, color=:black)
    plot!(fig, [0., Rneedle], [maximum(ztot)*1.15, maximum(ztot)*1.15], color=:black, lw=1)
    annotate!(fig, 0. + 0.25, maximum(ztot)*1.2, text(L"R_\mathrm{n}", :black))
    plot!(fig, [rtot[end], rtot[end]], [ztot[end], ztot[end]*1.2], color=:gray)
    plot!(fig, [rtot[1], rtot[1]], [ztot[end], ztot[end]*1.2], color=mypalette[4])

    # tangent
    plot!(fig, tanx, tany, color=:black, lw=1)
    dashed_arc!(fig, tanx[2], tany[2],r=0.15, θ1=0.0, θ2=pi/2.2, n=60, color=:black)
    annotate!(fig, tanx[2]+0.2, tany[2]+0.2, text(L"\psi", :black))

    # Rapex
    dashed_circle!(fig, 0.0, Rapex, Rapex; start_angle=1.3pi, end_angle=1.7pi, n=200, color=:black)
    radius_line!(fig, 0.0, Rapex, Rapex, 1.4pi; color=:black, lw=1)
    annotate!(fig, -0.5*Rapex+0.1, Rapex*0.5+0.05, text(L"R_\mathrm{apex}", :black))

    # material properties and pressure
    annotate!(fig, 0.45, 3.0, text(L"\rho_\mathrm{l}", :black))
    annotate!(fig, 1.5, 3.0, text(L"\rho_\mathrm{a}", :black))

    # g
    plot!(fig, [1.8, 1.8], [2.0, 1.0], arrow=:closed, lw=1, color=:black)
    annotate!(fig, 1.8+0.15, 1.5, text(L"\textbf{g}", :black))

    return fig
end


function plot_dropbounds_fn(dir_px; kwargs...)

    # --- JLD2 files ---
    jld2_files = filter(f -> endswith(f, ".jld2"),
                        readdir(dir_px; join=true))

    profiles_px = []

    for file in jld2_files
        @load file rpx Zpx
        push!(profiles_px, (rpx, Zpx, file))
    end

    # --- JPG image ---
    jpg_files = filter(f -> endswith(lowercase(f), ".jpg"),
                    readdir(dir_px; join=true))

    length(jpg_files) == 1 || error("Expected exactly one JPG")

    img = load(jpg_files[1])

    fig = plot(
          img,
          dpi=1000,      
          seriestype = :heatmap,
          aspect_ratio = :equal,
          xticks=false,
          yticks=false,
          size = (600, 600),
          framestyle = :none,
          left_margin = -10mm,
          bottom_margin = -10mm,
          right_margin = -10mm,
          top_margin = -10mm,
    )
 
    is = [12, 5, 6]
    for (i, (rs, Zs, name)) in enumerate(profiles_px)
        ii = is[i]
        xnew, ynew = reorder_dropletpx(rs, Zs)
        plot!(fig, xnew, ynew, label = false, lw=2., color=mypalette[ii])
    end

    return fig
end


function plot_lhood_fn(dir_output, dir_data; kwargs...)

    # load file
    file = joinpath(dir_output, "schematic.jld2")
    @load file r z s hits Rapex Rneedle

    isdir(dir_data) || error("Directory not found: $dir_data")

    # load undistorted
    filed = joinpath(dir_data, "dropletcoords_rhs_out.jld2")
    @load filed rr_sorted Zr_sorted

    rdf = rr_sorted
    zdf = Zr_sorted

    rdfull, zdfull, lastval = trim_trailing_repeats(rdf, zdf; keepmissing=false)

    Nrep = 16
    ds, sdfull = determine_arclength(rdfull, zdfull)
    rd, zd = trim_data(Nrep, rdfull, zdfull, sdfull)

    m = hcat([r, z]...)  

    idxs, dists = find_nearest_neighbors(m, rd, zd)
    figproj = plot(r, z, label=L"\mathrm{Model}", xlabel=L"r \, [\mathrm{mm}]", ylabel=L"z \, [\mathrm{mm}]", guidefontsize=16, tickfontsize=12, legendfontsize=12, color=mypalette[3], linewidth=2, size=(562.5,400), layout=grid(1,2, widths=(2.8/6,3.2/6)))
    scatter!(figproj, rd, zd, label=L"\mathrm{Data}", color=:black, markersize=2)
    # plot!(figproj, xlims=(-0.01, 2.5))
    plot!(figproj, xlim=(0, 2.5), ylim=(-4.1, 0.1), legend=:topright, legend_font_halign=:left)
    # Highlight closest neighbors
    xs_closest = r[idxs]
    ys_closest = z[idxs]
    for i in eachindex(rd)
        plot!(figproj, [rd[i], xs_closest[i]], [zd[i], ys_closest[i]], color=mypalette[4], label="")
    end
    plot!(figproj, [0.4, 1.1, 1.1, 0.4, 0.4], [0, 0, -1, -1, 0], color=:black, linestyle=:dash, label="", lw=1, subplot=1)
    # display(figproj)
    
    z_zoom_indices = findall(x -> x >= -1, z)
    r_zoom = r[z_zoom_indices]
    z_zoom = z[z_zoom_indices]
    plot!(figproj, r_zoom, z_zoom, xlim=(0.4,1.1), ylim=(-1.025, 0.025), xticks=(0.4:0.2:1), yticks=(0:-0.2:-1), legend=false, legend_font_halign=:left, color=mypalette[3], aspect_ratio=:equal, linewidth=2, subplot=2)
    # Find data points in zoom region
    rd_zoom_indices = findall(x -> x >= -1, zd)
    rd_zoom = rd[rd_zoom_indices]
    zd_zoom = zd[rd_zoom_indices]
    scatter!(figproj, rd_zoom, zd_zoom, xlabel=L"r \, [\mathrm{mm}]", ylabel=L"z \, [\mathrm{mm}]", color=:black, markersize=2, subplot=2, left_margin=10mm)
    for (i, orig_idx) in enumerate(rd_zoom_indices)
        # Use the original indices to get the correct nearest neighbor points
        if orig_idx <= length(idxs)
            plot!(figproj, [rd_zoom[i], xs_closest[orig_idx]], [zd_zoom[i], ys_closest[orig_idx]], color=mypalette[4], legend=false, subplot=2)
        end
    end
    if length(rd_zoom) >= 3
        # Use the third data point from the top
        orig_idx_third = rd_zoom_indices[3]
        if orig_idx_third <= length(idxs)
            annotate!(figproj, rd_zoom[3] + 0.12, zd_zoom[3], text(L"\delta_i", :black, 14), subplot=2)
        end
    end
    fig = figproj

    return fig
end


function plot_Rapex_fn(dir_data; kwargs...)

    # load file
    file = joinpath(dir_output, "schematic.jld2")
    @load file r z s hits Rapex Rneedle

    isdir(dir_data) || error("Directory not found: $dir_data")

    # load undistorted
    filed = joinpath(dir_data, "dropletcoords_out.jld2")
    @load filed rs_sorted Zs_sorted

    rptots, zptots = reorder_droplet(rs_sorted, Zs_sorted)
    rptotf, zptotf = trim_top_duplicates(rptots, zptots)

    Rapexs = 2*[0.4, 0.7, 1.0]
    fig = plot(xlabel=L"r \, [\mathrm{mm}]", ylabel=L"z \, [\mathrm{mm}]")
    plot!(fig, rptotf, zptotf, label=L"\mathrm{Experimental} \ \mathrm{data}", color=:black, lw=2)
    plot!(fig, xlims=(-1.5,4))
    for i = 1:3
        Rapex = Rapexs[i]
        mycircle!(fig, 0.0, minimum(zptotf)+Rapex, Rapex, start_angle=1.25π, end_angle=1.75π; color=mypalette[i+4], lw=2)#latexstring("R_\\mathrm{apex} = $(Rapex) \ \\mathrm{mm}"))
    end

   return fig

end