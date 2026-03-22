module TimCopPlotHelpers

using Graphs
using JLD2
using VortexDistributions

export DemoConfig,
    default_fixture_path,
    detect_fixture_demo,
    detect_validation_demo,
    ensure_example_plotting_modules,
    save_detection_figure,
    summarize_result

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
    @load sim_path sim
    @load sol_path sol

    psi = sol(t)
    X = sim.X
    x, y, z = X
    x = round.(x, digits = 3)
    y = round.(y, digits = 3)
    z = round.(z, digits = 3)
    Xr = (x, y, z)

    result = VortexDistributions.detect_vortices_3d(psi, x, y, z; n_itp = n_itp)
    connected = VortexDistributions.Detection3D.connect_vortex_ends(result, Xr)
    return psi, Xr, result, connected
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
