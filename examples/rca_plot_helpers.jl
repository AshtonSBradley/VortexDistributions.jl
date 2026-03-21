const RCA = VortexDistributions.RecursiveClusterAlgorithm

function classify_vortices(vortices, result)
    labels = fill("free vortex", length(vortices))
    colors = [vortex.qv > 0 ? "red" : "dodgerblue" for vortex in vortices]

    for dipole in result.dipoles
        for member in (dipole.vp, dipole.vn)
            idx = findfirst(isequal(member), vortices)
            if !isnothing(idx)
                labels[idx] = "dipole"
            end
        end
    end

    for (cluster_index, cluster) in pairs(result.positive_clusters)
        for member in cluster.vortices
            idx = findfirst(isequal(member), vortices)
            if !isnothing(idx)
                labels[idx] = "positive cluster $cluster_index"
            end
        end
    end

    for (cluster_index, cluster) in pairs(result.negative_clusters)
        for member in cluster.vortices
            idx = findfirst(isequal(member), vortices)
            if !isnothing(idx)
                labels[idx] = "negative cluster $cluster_index"
            end
        end
    end

    return labels, colors
end

function point_vortex_streamfunction(vortices, xgrid, ygrid; core = 0.03)
    ψ = zeros(length(ygrid), length(xgrid))
    core2 = core^2

    for (j, y) in pairs(ygrid)
        for (i, x) in pairs(xgrid)
            value = 0.0
            for vortex in vortices
                dx = x - vortex.xv
                dy = y - vortex.yv
                r2 = dx^2 + dy^2 + core2
                value -= vortex.qv * log(r2) / (4π)
            end
            ψ[j, i] = value
        end
    end

    return ψ
end

function plot_rca_connections!(plt, result)
    dipole_color = "mediumseagreen"
    positive_color = "red"
    negative_color = "dodgerblue"

    for dipole in result.dipoles
        plot!(
            plt,
            [dipole.vp.xv, dipole.vn.xv],
            [dipole.vp.yv, dipole.vn.yv],
            color = dipole_color,
            linewidth = 2.8,
            alpha = 0.9,
            label = "",
        )
    end

    for cluster in result.positive_clusters
        for edge in cluster.tree
            src = cluster.vortices[edge.src]
            dst = cluster.vortices[edge.dst]
            plot!(
                plt,
                [src.xv, dst.xv],
                [src.yv, dst.yv],
                color = positive_color,
                linewidth = 2.2,
                alpha = 0.9,
                label = "",
            )
        end
    end

    for cluster in result.negative_clusters
        for edge in cluster.tree
            src = cluster.vortices[edge.src]
            dst = cluster.vortices[edge.dst]
            plot!(
                plt,
                [src.xv, dst.xv],
                [src.yv, dst.yv],
                color = negative_color,
                linewidth = 2.2,
                alpha = 0.9,
                label = "",
            )
        end
    end
end

function plot_rca_example(vortices, result, labels, colors; Lx = 2.0, Ly = 2.0, title = "Recursive Cluster Algorithm Classification")
    gr()
    y_line_offset = 0.015 * Ly

    plt = plot(
        xlabel = "x",
        ylabel = "y",
        title = title,
        xlims = (-Lx / 2, Lx / 2),
        ylims = (-Ly / 2, Ly / 2),
        aspect_ratio = :equal,
        legend = :outertop,
        framestyle = :box,
        grid = true,
        background_color = RGB(0.99, 0.99, 0.97),
        foreground_color_grid = RGBA(0.7, 0.7, 0.7, 0.35),
        size = (900, 600),
    )

    xgrid = range(-Lx / 2, Lx / 2; length = 100)
    ygrid = range(-Ly / 2, Ly / 2; length = 100)
    ψ = point_vortex_streamfunction(vortices, xgrid, ygrid)

    contour!(
        plt,
        collect(xgrid),
        collect(ygrid),
        ψ,
        levels = 40,
        c = :thermal,
        linewidth = 1.8,
        alpha = 0.65,
        colorbar = false,
        label = "",
    )

    plot_rca_connections!(plt, result)

    for (idx, vortex) in pairs(vortices)
        annotate!(
            plt,
            vortex.xv,
            vortex.yv + (vortex.qv > 0 ? y_line_offset / 2 : 7y_line_offset / 10),
            text(vortex.qv > 0 ? "+" : "-", vortex.qv > 0 ? 15 : 30, colors[idx], :center, "Computer Modern"),
        )
    end

    legend_entries = [
        ("positive cluster", "red"),
        ("negative cluster", "dodgerblue"),
        ("vortex dipole", "mediumseagreen"),
    ]

    for (label, color) in legend_entries
        if label == "positive cluster" && isempty(result.positive_clusters)
            continue
        elseif label == "negative cluster" && isempty(result.negative_clusters)
            continue
        elseif label == "vortex dipole" && isempty(result.dipoles)
            continue
        end
        plot!(
            plt,
            [NaN, NaN],
            [NaN, NaN],
            color = color,
            linewidth = 4.5,
            label = label,
        )
    end

    return plt
end

function print_rca_summary(result, plot_path)
    println("Dipoles: ", length(result.dipoles))
    println("Positive clusters: ", length(result.positive_clusters))
    println("Negative clusters: ", length(result.negative_clusters))
    println("Free vortices: ", length(result.free_vortices))
    println("Saved plot: ", plot_path)

    for (i, cluster) in pairs(result.positive_clusters)
        println("Positive cluster $i has $(length(cluster.vortices)) vortices")
    end

    for (i, cluster) in pairs(result.negative_clusters)
        println("Negative cluster $i has $(length(cluster.vortices)) vortices")
    end
end
