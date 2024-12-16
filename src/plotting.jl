function plot_unconnected_vorts_mat(vort_lines, vort_loops, vort_rings; linewidth=5)
    colors = distinguishable_colors(length(vort_lines) + length(vort_loops) + length(vort_rings));

    all_vorts = []
    all_colors = []

    for i in eachindex(vort_lines)
        for j in 1:size(vort_lines[i])[1]
            v_sort = vort_lines[i][j, :]
            push!(all_vorts, (v_sort[1], v_sort[2], v_sort[3]))
            push!(all_colors, colors[i])
        end
        push!(all_vorts, (NaN, NaN, NaN))
        push!(all_colors, colors[i])
    end

    # for i in eachindex(vort_loops)
    #     for j in eachindex(vort_loops[i])
    #         v_sort = vort_loops[i][j]
    #         push!(all_vorts, (v_sort[1], v_sort[2], v_sort[3]))
    #         push!(all_colors, colors[i + length(vort_lines)])
    #     end
    #     push!(all_vorts, (NaN, NaN, NaN))
    #     push!(all_colors, colors[i + length(vort_lines)])
    # end

    for i in eachindex(vort_loops)
        for j in 1:size(vort_loops[i])[1]
            v_sort = vort_loops[i][j, :]
            push!(all_vorts, (v_sort[1], v_sort[2], v_sort[3]))
            push!(all_colors, colors[i + length(vort_lines)])
        end
        push!(all_vorts, (NaN, NaN, NaN))
        push!(all_colors, colors[i + length(vort_lines)])
    end

    # for i in eachindex(vort_rings)
    #     for j in eachindex(vort_rings[i])
    #         v_sort = vort_rings[i][j]
    #         push!(all_vorts, (v_sort[1], v_sort[2], v_sort[3]))
    #         push!(all_colors, colors[i + length(vort_lines) + length(vort_loops)])
    #     end
    #     push!(all_vorts, (NaN, NaN, NaN))
    #     push!(all_colors, colors[i + length(vort_lines) + length(vort_loops)])
    # end

    for i in eachindex(vort_rings)
        for j in 1:size(vort_rings[i])[1]
            v_sort = vort_rings[i][j, :]
            push!(all_vorts, (v_sort[1], v_sort[2], v_sort[3]))
            push!(all_colors, colors[i + length(vort_lines) + length(vort_loops)])
        end
        push!(all_vorts, (NaN, NaN, NaN))
        push!(all_colors, colors[i + length(vort_lines) + length(vort_loops)])
    end
    
    lines!(all_vorts, color=all_colors, linewidth=linewidth)
end