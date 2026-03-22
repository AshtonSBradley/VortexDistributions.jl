module TimCopPlotHelpers

using FLoops
using Graphs
using JLD2
using VortexDistributions

include(joinpath(@__DIR__, "timcop_ccma.jl"))
using .TimCopCCMA: CCMA, filter

export DemoConfig,
    coords_to_path_matrices,
    default_fixture_path,
    detect_fixture_demo,
    detect_validation_demo,
    ensure_example_plotting_modules,
    graphplot_detection!,
    load_validation_state,
    plot_iso,
    plot_unconnected_vorts_mat!,
    save_detection_figure,
    summarize_result,
    vortex_ccma_j

Base.@kwdef struct DemoConfig
    n_itp::Int = 1
    output_path::String = joinpath(@__DIR__, "timcop_detection_demo.png")
    isovalue::Float64 = 0.35
    isorange::Float64 = 0.05
end

default_fixture_path() = normpath(joinpath(@__DIR__, "..", "test", "3d", "box_vorts.jld2"))

function ensure_example_plotting_modules()
    missing = String[]
    for pkg in ("GLMakie", "GraphMakie")
        if Base.find_package(pkg) === nothing
            push!(missing, pkg)
        end
    end

    if !isempty(missing)
        quoted = join(["\"$pkg\"" for pkg in missing], ", ")
        error(
            "This usage example requires $(join(missing, " and ")) in the active environment.\n" *
            "Install them with:\n" *
            "  julia --project=. -e 'using Pkg; Pkg.add([$quoted])'"
        )
    end

    Core.eval(@__MODULE__, :(using GLMakie))
    Core.eval(@__MODULE__, :(using GraphMakie))
    glmakie = Base.invokelatest(getfield, @__MODULE__, :GLMakie)
    graphmakie = Base.invokelatest(getfield, @__MODULE__, :GraphMakie)
    return glmakie, graphmakie
end

function detect_fixture_demo(; n_itp::Int = 1)
    fixture_path = default_fixture_path()
    @load fixture_path psi_tubes1 X
    psi = psi_tubes1
    x, y, z = X
    result = VortexDistributions.detect_vortices_3d(psi, x, y, z; n_itp = n_itp)
    connected = VortexDistributions.Detection3D.connect_vortex_ends(result, X)
    return psi, X, result, connected
end

function detect_validation_demo(sim_path::AbstractString, sol_path::AbstractString; t::Int = 20, n_itp::Int = 4)
    psi, Xr = load_validation_state(sim_path, sol_path; t = t)
    x, y, z = Xr
    result = VortexDistributions.detect_vortices_3d(psi, x, y, z; n_itp = n_itp)
    connected = VortexDistributions.Detection3D.connect_vortex_ends(result, Xr)
    return psi, Xr, result, connected
end

function load_validation_state(sim_path::AbstractString, sol_path::AbstractString; t::Int = 20)
    @load sim_path sim
    @load sol_path sol

    psi = sol(t)
    X = sim.X
    x, y, z = X
    x = round.(x, digits = 3)
    y = round.(y, digits = 3)
    z = round.(z, digits = 3)
    Xr = (x, y, z)
    return psi, Xr
end

function summarize_result(io::IO, result, connected; label::AbstractString)
    println(io, label)
    println(io, "  graph vertices: $(length(result.coords))")
    println(io, "  graph edges: $(Graphs.ne(result.graph))")
    println(io, "  vortex lines: $(length(result.vortex_lines))")
    println(io, "  vortex loops: $(length(result.vortex_loops))")
    println(io, "  vortex rings: $(length(result.vortex_rings))")
    println(io, "  connected structures: $(length(connected))")
    println(io, "  periodic edges: $(length(result.periodic_edges))")
end

function _line_segments(result)
    return (
        lines = result.vortex_lines,
        loops = result.vortex_loops,
        rings = result.vortex_rings,
    )
end

function _plot_structure!(backend, scene, coords, segments; color, linewidth)
    for segment in segments
        points = coords[segment]
        xs = first.(points)
        ys = getindex.(points, 2)
        zs = getindex.(points, 3)
        backend.lines!(scene, xs, ys, zs; color = color, linewidth = linewidth)
    end
end

function plot_iso(psi, X; visible = true, isovalue = 0.7, isorange = 0.1, show_axis = true)
    mk, _ = ensure_example_plotting_modules()
    density = abs2.(psi)
    density ./= maximum(density)
    x, y, z = X

    fig = mk.Figure(size = (1200, 900))
    scene = mk.LScene(fig[1, 1], show_axis = show_axis)
    mk.volume!(
        scene,
        (x[1], x[end]),
        (y[1], y[end]),
        (z[1], z[end]),
        density;
        algorithm = :iso,
        isovalue = isovalue,
        isorange = isorange,
        transparency = true,
        visible = visible,
        colormap = :plasma,
    )

    return fig, scene
end

function graphplot_detection!(scene, graph, coords; edge_width = 1, node_size = 2)
    _, gm = ensure_example_plotting_modules()
    gm.graphplot!(scene, graph; layout = (_) -> coords, edge_width = edge_width, node_size = node_size, show_axis = false)
    return scene
end

function coords_to_path_matrices(coords, paths)
    return [reduce(vcat, transpose.(coords[path])) for path in paths]
end

function vortex_ccma_j(vort_lines_mat, vort_loops_mat, vort_rings_mat; distrib = "hanning", w_ma = 2, w_cc = 3,
                       itr = 1, cc_mode = true)
    ccma_inst = CCMA(w_ma = w_ma, w_cc = w_cc, distrib = distrib)
    vort_lines_ccma = Vector{Any}(undef, length(vort_lines_mat))

    @floop for (i, v_l) in enumerate(vort_lines_mat)
        try
            vort_lines_ccma[i] = filter(ccma_inst, v_l, mode = "fill_boundary", cc_mode = cc_mode)
        catch
            vort_lines_ccma[i] = v_l
        end
    end

    vort_loops_ccma = Vector{Any}(undef, length(vort_loops_mat))
    @floop for (i, v_l) in enumerate(vort_loops_mat)
        try
            vort_loops_ccma[i] = filter(ccma_inst, v_l, mode = "fill_boundary", cc_mode = cc_mode)
        catch
            vort_loops_ccma[i] = v_l
        end
    end

    vort_rings_ccma = Vector{Any}(undef, length(vort_rings_mat))
    @floop for (i, v_l) in enumerate(vort_rings_mat)
        try
            vort_rings_ccma[i] = filter(ccma_inst, v_l, mode = "wrapping", cc_mode = cc_mode)
        catch
            vort_rings_ccma[i] = v_l
        end
    end

    return vort_lines_ccma, vort_loops_ccma, vort_rings_ccma
end

function plot_unconnected_vorts_mat!(scene, vort_lines, vort_loops, vort_rings; linewidth = 5, mono = false,
                                     color = :black, repeat_ring_ends = false)
    mk, _ = ensure_example_plotting_modules()
    colors = mk.distinguishable_colors(length(vort_lines) + length(vort_loops) + length(vort_rings))

    all_vorts = NTuple{3, Float64}[]
    all_colors = Any[]

    for i in eachindex(vort_lines)
        for j in axes(vort_lines[i], 1)
            v_sort = vort_lines[i][j, :]
            push!(all_vorts, (v_sort[1], v_sort[2], v_sort[3]))
            push!(all_colors, colors[i])
        end
        push!(all_vorts, (NaN, NaN, NaN))
        push!(all_colors, colors[i])
    end

    for i in eachindex(vort_loops)
        for j in axes(vort_loops[i], 1)
            v_sort = vort_loops[i][j, :]
            push!(all_vorts, (v_sort[1], v_sort[2], v_sort[3]))
            push!(all_colors, colors[i + length(vort_lines)])
        end
        push!(all_vorts, (NaN, NaN, NaN))
        push!(all_colors, colors[i + length(vort_lines)])
    end

    for i in eachindex(vort_rings)
        for j in axes(vort_rings[i], 1)
            v_sort = vort_rings[i][j, :]
            push!(all_vorts, (v_sort[1], v_sort[2], v_sort[3]))
            push!(all_colors, colors[i + length(vort_lines) + length(vort_loops)])
        end
        if repeat_ring_ends && size(vort_rings[i], 1) > 0
            push!(all_vorts, (vort_rings[i][1, 1], vort_rings[i][1, 2], vort_rings[i][1, 3]))
            push!(all_colors, colors[i + length(vort_lines) + length(vort_loops)])
        end
        push!(all_vorts, (NaN, NaN, NaN))
        push!(all_colors, colors[i + length(vort_lines) + length(vort_loops)])
    end

    line_colors = mono ? color : all_colors
    mk.lines!(scene, all_vorts; color = line_colors, linewidth = linewidth)
    return scene
end

function _save_detection_figure_loaded(mk, gm, psi, X, result, output_path, isovalue, isorange)
    mkpath(dirname(output_path))

    x, y, z = X
    density = abs2.(psi)
    density ./= maximum(density)

    fig = mk.Figure(size = (1200, 900))
    mk.Label(fig[0, 1], "Detected vortex structures", fontsize = 22)
    mk.Label(fig[0, 2], "Detection graph (GraphMakie)", fontsize = 22)

    filament_scene = mk.LScene(fig[1, 1], show_axis = true)
    graph_scene = mk.LScene(fig[1, 2], show_axis = true)

    mk.volume!(
        filament_scene,
        (x[1], x[end]),
        (y[1], y[end]),
        (z[1], z[end]),
        density;
        algorithm = :iso,
        isovalue = isovalue,
        isorange = isorange,
        transparency = true,
        colormap = :viridis,
    )

    structures = _line_segments(result)
    _plot_structure!(mk, filament_scene, result.coords, structures.lines; color = :dodgerblue3, linewidth = 4)
    _plot_structure!(mk, filament_scene, result.coords, structures.loops; color = :darkorange2, linewidth = 4)
    _plot_structure!(mk, filament_scene, result.coords, structures.rings; color = :crimson, linewidth = 4)

    gm.graphplot!(
        graph_scene,
        result.graph;
        layout = (_) -> result.coords,
        edge_width = 1.5,
        node_size = 5,
    )

    mk.save(output_path, fig)
    return output_path
end

function save_detection_figure(
    psi,
    X,
    result,
    output_path::AbstractString;
    isovalue::Float64 = 0.35,
    isorange::Float64 = 0.05,
)
    mk, gm = ensure_example_plotting_modules()
    return Base.invokelatest(_save_detection_figure_loaded, mk, gm, psi, X, result, output_path, isovalue, isorange)
end

end
