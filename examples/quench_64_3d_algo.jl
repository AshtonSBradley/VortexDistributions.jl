using JLD2
using Graphs
using VortexDistributions

include(joinpath(@__DIR__, "timcop_plot_helpers.jl"))
using .TimCopPlotHelpers

"""
Upstream-style 64^3 3D detection and plotting example adapted from
`timcop/VortexDetection.jl/src/examples/quench_64_3d_algo.jl`.

Usage:

    julia --project=. examples/quench_64_3d_algo.jl sim_64.jld2 sol_64.jld2 [output_png] [t] [n_itp] [n_itr]

Defaults:
- `output_png = examples/quench_64_3d_algo.png`
- `t = 20`
- `n_itp = 4`
- `n_itr = 10`

This script intentionally requires `GLMakie` and `GraphMakie`, but keeps them
out of the main package dependency graph.
"""

function parse_args(args)
    if length(args) < 2
        error("Usage: julia --project=. examples/quench_64_3d_algo.jl sim_64.jld2 sol_64.jld2 [output_png] [t] [n_itp] [n_itr]")
    end

    sim_path = args[1]
    sol_path = args[2]
    output_path = length(args) >= 3 ? args[3] : joinpath(@__DIR__, "quench_64_3d_algo.png")
    t = length(args) >= 4 ? parse(Int, args[4]) : 20
    n_itp = length(args) >= 5 ? parse(Int, args[5]) : 4
    n_itr = length(args) >= 6 ? parse(Int, args[6]) : 10

    return sim_path, sol_path, output_path, t, n_itp, n_itr
end

function main(args)
    sim_path, sol_path, output_path, t, n_interpolate, n_itr = parse_args(args)

    psi, X = load_validation_state(sim_path, sol_path; t = t)
    x, y, z = X

    fig, lscene = plot_iso(
        psi,
        X;
        show_axis = true,
        isovalue = 0.5,
        isorange = 0.15,
        visible = true,
    )

    println("Running upstream-style TimCop 3D workflow")
    println("  time index: $t")
    println("  interpolation factor: $n_interpolate")
    println("  smoothing iterations: $n_itr")

    result = @timed VortexDistributions.full_algorithm(psi, x, y, z; n_itp = n_interpolate)
    g, vort_lines, vort_loops, vort_rings, vorts_coords, periodic_edges = result.value
    println("  detection time: $(round(result.time; digits = 3)) s")

    graphplot_detection!(lscene, g, vorts_coords; edge_width = 1, node_size = 2)

    vort_lines_coords = coords_to_path_matrices(vorts_coords, vort_lines)
    vort_loops_coords = coords_to_path_matrices(vorts_coords, vort_loops)
    vort_rings_coords = coords_to_path_matrices(vorts_coords, vort_rings)

    smoothing = @timed vortex_ccma_j(
        vort_lines_coords,
        vort_loops_coords,
        vort_rings_coords;
        w_ma = 3,
        w_cc = 2,
        distrib = "hanning",
    )
    vort_lines_ccma_j, vort_loops_ccma_j, vort_rings_ccma_j = smoothing.value
    println("  initial smoothing pass: $(round(smoothing.time; digits = 3)) s")

    for _ in 1:n_itr
        vort_lines_ccma_j, vort_loops_ccma_j, vort_rings_ccma_j = vortex_ccma_j(
            vort_lines_ccma_j,
            vort_loops_ccma_j,
            vort_rings_ccma_j;
            w_ma = 3,
            w_cc = 2,
            distrib = "hanning",
        )
    end

    plot_unconnected_vorts_mat!(
        lscene,
        vort_lines_ccma_j,
        vort_loops_ccma_j,
        vort_rings_ccma_j;
        linewidth = 3,
        mono = true,
        repeat_ring_ends = true,
        color = :blue,
    )

    connected = VortexDistributions.Detection3D.connect_vortex_ends(vorts_coords, vort_lines, vort_loops, vort_rings, X)
    println("Quench 64 3D detection summary")
    println("  graph vertices: $(length(vorts_coords))")
    println("  graph edges: $(Graphs.ne(g))")
    println("  vortex lines: $(length(vort_lines))")
    println("  vortex loops: $(length(vort_loops))")
    println("  vortex rings: $(length(vort_rings))")
    println("  connected structures: $(length(connected))")
    println("  periodic edges: $(length(periodic_edges))")

    mk, _ = ensure_example_plotting_modules()
    mk.save(output_path, fig)
    println("  wrote plot: $output_path")
end

main(ARGS)
