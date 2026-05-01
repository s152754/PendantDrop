function llhooddist_fig(rm, zm, rd, zd)

        # rm, zm = modelfn(log.([72.0, 0.7]))  # example parameters: sigma, R
        m = hcat([rm, zm]...)  

        idxs, dists = find_nearest_neighbors(m, rd, zd)
        figproj = plot(rm, zm, label=L"\mathrm{Model}", color=mypalette[3], aspect_ratio=:equal, lw=3)
        scatter!(figproj, rd, zd, label=L"\mathrm{Data}", color=:black, markersize=3)
        # plot!(figproj, xlims=(-0.01, 2.5))
        plot!(figproj, xlims=(-0.01, 2.5), ylims=(-4.0, 0.02), size=(415,620), legend=:topright, legend_font_halign=:left)
        # Highlight closest neighbors
        xs_closest = rm[idxs]
        ys_closest = zm[idxs]
        for i in eachindex(rd)
            plot!(figproj, [rd[i], xs_closest[i]], [zd[i], ys_closest[i]], color=mypalette[4], label="")
        end
        display(figproj)

        # Draw rectangle for zoom region
        zoom_xlims = (0.4, 1.1)
        zoom_ylims = (-1., 0.02)
        plot!(figproj, [zoom_xlims[1], zoom_xlims[2], zoom_xlims[2], zoom_xlims[1], zoom_xlims[1]],
            [zoom_ylims[1], zoom_ylims[1], zoom_ylims[2], zoom_ylims[2], zoom_ylims[1]],
            color=:black, linestyle=:dash, label="", lw=2)

        # Create inset plot
        
        inset_plot = plot(xlims=zoom_xlims, ylims=zoom_ylims, aspect_ratio=:equal, framestyle=:none, size=())
        scatter!(inset_plot, rd, zd, color=:black, markersize=3, label="")
        plot!(inset_plot, size=(415, 500))
        plot!(inset_plot, rm, zm, color=mypalette[3], label="", lw=3)
        plot!(inset_plot, xlims=zoom_xlims, ylims=zoom_ylims, aspect_ratio=:equal,
            grid=false, yticks=false, xticks=false)
        plot!(inset_plot, [zoom_xlims[1], zoom_xlims[2], zoom_xlims[2], zoom_xlims[1], zoom_xlims[1]],
            [zoom_ylims[1], zoom_ylims[1], zoom_ylims[2], zoom_ylims[2], zoom_ylims[1]],
            color=:black, linestyle=:dash, label="", lw=3)


        # Draw distance vectors and annotate

        for i in eachindex(rd)
            global hits
            if zoom_xlims[1] ≤ rd[i] ≤ zoom_xlims[2] && zoom_ylims[1] ≤ zd[i] ≤ zoom_ylims[2]
                dx = rm[idxs[i]] - rd[i]
                dy = zm[idxs[i]] - zd[i]
                # quiver!(inset_plot, [rd[i]], [zd[i]], quiver=([dx], [dy]), color=:black, arrow=:closed, quiver_scale=0.8)
                plot!(inset_plot, [rd[i], rd[i] + dx], [zd[i], zd[i] + dy], color=mypalette[4], linewidth=1, label="")      # smaller lw => smaller-looking
                midx = rd[i] + dx/2
                midy = zd[i] + dy/2
                
                # Normal vector
                nx = -dy / hypot(dx, dy)
                ny = dx / hypot(dx, dy)

                offset = 0.05  # adjust as needed
                normalx = midx + offset * nx
                normaly = midy + offset * ny
                hits += 1

                # annotate!(inset_plot, midx, midy, text(L"\delta_i", :black, 10))
                if hits == 1
                    annotate!(inset_plot, normalx, normaly, text(L"\delta_i", :black, 14))
                end
                # break
            end
        end

        # Combine main plot and inset
        final_plot = plot(figproj, inset_plot, layout=@layout([a{0.5w} b{0.4w}]))
        # plot!(final_plot, guidefont=18, tickfont=14, legendfont=14, bottom_margin=4mm, left_margin=-10mm)
        plot!(final_plot, size=(600, 415), labelfontsize=20, tickfontsize=16, annotationfontsize=17, legendfontsize=17)

        return final_plot
        
end
##
function llhooddist_figgg(rm, zm, rd, zd)

        # rm, zm = modelfn(log.([72.0, 0.7]))  # example parameters: sigma, R
        m = hcat([rm, zm]...)  

        idxs, dists = find_nearest_neighbors(m, rd, zd)
        figproj = plot(rm, zm, label=L"\mathrm{Model}", xlabel=L"r \, [\mathrm{mm}]", ylabel=L"z \, [\mathrm{mm}]", guidefontsize=16, tickfontsize=12, legendfontsize=12, color=mypalette[3], linewidth=2, size=(562.5,400), layout=grid(1,2, widths=(2.8/6,3.2/6)))
        scatter!(figproj, rd, zd, label=L"\mathrm{Data}", color=:black, markersize=2)
        # plot!(figproj, xlims=(-0.01, 2.5))
        plot!(figproj, xlim=(0, 2.5), ylim=(-4.1, 0.1), legend=:topright, legend_font_halign=:left)
        # Highlight closest neighbors
        xs_closest = rm[idxs]
        ys_closest = zm[idxs]
        for i in eachindex(rd)
            plot!(figproj, [rd[i], xs_closest[i]], [zd[i], ys_closest[i]], color=mypalette[4], label="")
        end
        plot!(figproj, [0.4, 1.1, 1.1, 0.4, 0.4], [0, 0, -1, -1, 0], color=:black, linestyle=:dash, label="", lw=1, subplot=1)
        # display(figproj)
        
        zm_zoom = zm[zm .>= -1]
        index_zm = findall(x -> x == zm_zoom[1], zm)
        rm_zoom = rm[index_zm[1]:end]
        zd_zoom = zd[zd .>= -1]
        plot!(figproj, rm_zoom, zm_zoom, xlim=(0.4,1.1), ylim=(-1.025, 0.025), xticks=(0.4:0.2:1), yticks=(0:-0.2:-1), legend=false, legend_font_halign=:left, color=mypalette[3], aspect_ratio=:equal, linewidth=2, subplot=2)
        scatter!(figproj, rd[end-length(zd_zoom)+1:end], zd_zoom, xlabel=L"r \, [\mathrm{mm}]", ylabel=L"z \, [\mathrm{mm}]", color=:black, markersize=2, subplot=2, left_margin=10mm)
        for i in eachindex(zd_zoom)
            plot!(figproj, [rd[end-i+1], xs_closest[end-i+1]], [zd[end-i+1], ys_closest[end-i+1]], color=mypalette[4], legend=false, subplot=2)
        end
        annotate!(figproj, mean((rd[end-(length(zd_zoom)-1)], xs_closest[end-(length(zd_zoom)-1)])), zd[end-(length(zd_zoom)-1)]+(zd[end-(length(zd_zoom)-1)]-ys_closest[end-(length(zd_zoom)-1)]), text(L"\delta_i", :black, 14), subplot=2)
        final_plot = figproj
        # inset_plot = lens([0.4, 1.1], [-1,0])

        # # Draw rectangle for zoom region
        # zoom_xlims = (0.4, 1.1)
        # zoom_ylims = (-1., 0.02)
        # plot!(figproj, [zoom_xlims[1], zoom_xlims[2], zoom_xlims[2], zoom_xlims[1], zoom_xlims[1]],
        #     [zoom_ylims[1], zoom_ylims[1], zoom_ylims[2], zoom_ylims[2], zoom_ylims[1]],
        #     color=:black, linestyle=:dash, label="", lw=2)

        # # Create inset plot
        
        # inset_plot = plot(xlims=zoom_xlims, ylims=zoom_ylims, aspect_ratio=:equal, framestyle=:none, size=())
        # scatter!(inset_plot, rd, zd, color=:black, markersize=3, label="")
        # plot!(inset_plot, size=(415, 500))
        # plot!(inset_plot, rm, zm, color=mypalette[3], label="", lw=3)
        # plot!(inset_plot, xlims=zoom_xlims, ylims=zoom_ylims, aspect_ratio=:equal,
        #     grid=false, yticks=false, xticks=false)
        # plot!(inset_plot, [zoom_xlims[1], zoom_xlims[2], zoom_xlims[2], zoom_xlims[1], zoom_xlims[1]],
        #     [zoom_ylims[1], zoom_ylims[1], zoom_ylims[2], zoom_ylims[2], zoom_ylims[1]],
        #     color=:black, linestyle=:dash, label="", lw=3)


        # # Draw distance vectors and annotate

        # for i in eachindex(rd)
        #     global hits
        #     if zoom_xlims[1] ≤ rd[i] ≤ zoom_xlims[2] && zoom_ylims[1] ≤ zd[i] ≤ zoom_ylims[2]
        #         dx = rm[idxs[i]] - rd[i]
        #         dy = zm[idxs[i]] - zd[i]
        #         # quiver!(inset_plot, [rd[i]], [zd[i]], quiver=([dx], [dy]), color=:black, arrow=:closed, quiver_scale=0.8)
        #         plot!(inset_plot, [rd[i], rd[i] + dx], [zd[i], zd[i] + dy], color=mypalette[4], linewidth=1, label="")      # smaller lw => smaller-looking
        #         midx = rd[i] + dx/2
        #         midy = zd[i] + dy/2
                
        #         # Normal vector
        #         nx = -dy / hypot(dx, dy)
        #         ny = dx / hypot(dx, dy)

        #         offset = 0.05  # adjust as needed
        #         normalx = midx + offset * nx
        #         normaly = midy + offset * ny
        #         hits += 1

        #         # annotate!(inset_plot, midx, midy, text(L"\delta_i", :black, 10))
        #         if hits == 1
        #             annotate!(inset_plot, normalx, normaly, text(L"\delta_i", :black, 14))
        #         end
        #         # break
        #     end
        # end

        # # Combine main plot and inset
        # final_plot = figproj
        # final_plot2 = plot(figproj, inset_plot, layout=@layout([a{0.5w} b{0.4w}]))
        # # plot!(final_plot, guidefont=18, tickfont=14, legendfont=14, bottom_margin=4mm, left_margin=-10mm)
        # plot!(final_plot, size=(600, 415), labelfontsize=20, tickfontsize=16, annotationfontsize=17, legendfontsize=17)

        return final_plot, zm_zoom, index_zm, zd_zoom#, rm_zoom
        
end
##

function priorRapex_fig(Rapexs, rptotf, zptotf)

        # rm, zm = modelfn(log.([72.0, 0.7]))  # example parameters: sigma, R
        m = hcat([rm, zm]...)  

        idxs, dists = find_nearest_neighbors(m, rd, zd)
        plot!(figproj, rm, zm, label=L"\mathrm{model}", color=mypalette[3], aspect_ratio=:equal, legend=:topright, lw=3)
        scatter!(figproj, rd, zd, label=L"\mathrm{data}", color=:black, markersize=4)
        plot!(figproj, xlims=(-0.01, 2.5))
        # Highlight closest neighbors
        xs_closest = rm[idxs]
        ys_closest = zm[idxs]
        for i in eachindex(rd)
            plot!(figproj, [rd[i], xs_closest[i]], [zd[i], ys_closest[i]], color=mypalette[4], label="")
        end
        display(figproj)


        # Draw rectangle for zoom region
        zoom_xlims = (0.4, 1.1)
        zoom_ylims = (-1., 0.02)
        plot!(figproj, [zoom_xlims[1], zoom_xlims[2], zoom_xlims[2], zoom_xlims[1], zoom_xlims[1]],
            [zoom_ylims[1], zoom_ylims[1], zoom_ylims[2], zoom_ylims[2], zoom_ylims[1]],
            color=:black, linestyle=:dash, label="", lw=2)

        # Create inset plot
        inset_plot = plot(xlims=zoom_xlims, ylims=zoom_ylims, aspect_ratio=:equal, framestyle=:none)
        scatter!(inset_plot, rd, zd, color=:black, markersize=4, label="")
        plot!(inset_plot, rm, zm, color=mypalette[3], label="", lw=3)
        plot!(inset_plot, xlims=zoom_xlims, ylims=zoom_ylims, aspect_ratio=:equal,
            grid=false, yticks=false, xticks=false)
        plot!(inset_plot, [zoom_xlims[1], zoom_xlims[2], zoom_xlims[2], zoom_xlims[1], zoom_xlims[1]],
            [zoom_ylims[1], zoom_ylims[1], zoom_ylims[2], zoom_ylims[2], zoom_ylims[1]],
            color=:black, linestyle=:dash, label="", lw=3)


        # Draw distance vectors and annotate

        for i in eachindex(rd)
            global hits
            if zoom_xlims[1] ≤ rd[i] ≤ zoom_xlims[2] && zoom_ylims[1] ≤ zd[i] ≤ zoom_ylims[2]
                dx = rm[idxs[i]] - rd[i]
                dy = zm[idxs[i]] - zd[i]
                # quiver!(inset_plot, [rd[i]], [zd[i]], quiver=([dx], [dy]), color=:black, arrow=:closed, quiver_scale=0.8)
                plot!(inset_plot, [rd[i], rd[i] + dx], [zd[i], zd[i] + dy], color=mypalette[4], linewidth=1, label="")      # smaller lw => smaller-looking
                midx = rd[i] + dx/2
                midy = zd[i] + dy/2
                
                # Normal vector
                nx = -dy / hypot(dx, dy)
                ny = dx / hypot(dx, dy)

                offset = 0.05  # adjust as needed
                normalx = midx + offset * nx
                normaly = midy + offset * ny
                hits += 1

                # annotate!(inset_plot, midx, midy, text(L"\delta_i", :black, 10))
                if hits == 1
                    annotate!(inset_plot, normalx, normaly, text(L"\delta_i", :black, 14))
                end
                # break
            end
        end

        # Combine main plot and inset
        final_plot = plot(figproj, inset_plot, layout=@layout([a{0.5w} b{0.4w}]))
        plot!(final_plot, guidefont=18, tickfont=14, legendfont=14, bottom_margin=4mm, left_margin=-10mm)

        return final_plot
        
end