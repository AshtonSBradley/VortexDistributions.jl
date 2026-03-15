using Graphs
using JLD2
using VortexDistributions

"""
Run a deeper 3D validation pass against simulation output generated from the
local `examples/generate_3d_validation_data.jl` workflow.

Usage:

    julia --project=. examples/timcop_deep_test.jl path/to/sim_64.jld2 path/to/sol_64.jld2 [t] [n_itp]

Arguments:
- `sim_64.jld2`: expected to contain `sim.X`
- `sol_64.jld2`: expected to contain `sol`, a callable solution object
- `t`: optional simulation time index, default `20`
- `n_itp`: optional interpolation factor, default `4`
"""

function main(args)
    if length(args) < 2
        error("Usage: julia --project=. examples/timcop_deep_test.jl sim_64.jld2 sol_64.jld2 [t] [n_itp]")
    end

    sim_path = args[1]
    sol_path = args[2]
    t = length(args) >= 3 ? parse(Int, args[3]) : 20
    n_itp = length(args) >= 4 ? parse(Int, args[4]) : 4

    @load sim_path sim
    @load sol_path sol

    psi = sol(t)
    X = sim.X
    x, y, z = X
    x = round.(x, digits = 3)
    y = round.(y, digits = 3)
    z = round.(z, digits = 3)

    result = detect_vortices_3d(psi, x, y, z; n_itp = n_itp)
    connected = VortexDistributions.Detection3D.connect_vortex_ends(result, X)

    println("3D detection summary")
    println("  time index: $t")
    println("  interpolation factor: $n_itp")
    println("  graph vertices: $(Graphs.nv(result.graph))")
    println("  graph edges: $(Graphs.ne(result.graph))")
    println("  filament lines: $(length(result.vortex_lines))")
    println("  filament loops: $(length(result.vortex_loops))")
    println("  filament rings: $(length(result.vortex_rings))")
    println("  connected vortex structures: $(length(connected))")
    println("  periodic edges: $(length(result.periodic_edges))")

    output_path = joinpath(dirname(sol_path), "vortex_detection_result_t$(t)_n$(n_itp).jld2")
    @save output_path result connected
    println("  wrote result bundle: $output_path")
end

main(ARGS)
