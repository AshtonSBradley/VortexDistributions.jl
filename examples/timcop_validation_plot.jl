include(joinpath(@__DIR__, "timcop_plot_helpers.jl"))
using .TimCopPlotHelpers

"""
Load the generated 64^3 validation data, run the integrated Tim Cop detector,
and save a GLMakie/GraphMakie visualization of the detected vortex structures.

Usage:

    julia --project=. examples/timcop_validation_plot.jl sim_64.jld2 sol_64.jld2 [output_png] [t] [n_itp]

Defaults:
- `output_png = examples/timcop_validation_plot.png`
- `t = 20`
- `n_itp = 4`

This complements `examples/generate_3d_validation_data.jl` and
`examples/timcop_deep_test.jl`, but remains outside the package test suite.
"""

function parse_args(args)
    positional = String[]
    for arg in args
        push!(positional, arg)
    end

    if length(positional) < 2
        error("Usage: julia --project=. examples/timcop_validation_plot.jl sim_64.jld2 sol_64.jld2 [output_png] [t] [n_itp]")
    end

    sim_path = positional[1]
    sol_path = positional[2]
    output_path = length(positional) >= 3 ? positional[3] : joinpath(@__DIR__, "timcop_validation_plot.png")
    t = length(positional) >= 4 ? parse(Int, positional[4]) : 20
    n_itp = length(positional) >= 5 ? parse(Int, positional[5]) : 4

    return sim_path, sol_path, DemoConfig(; output_path, n_itp), t
end

function main(args)
    sim_path, sol_path, config, t = parse_args(args)
    psi, X, result, connected = detect_validation_demo(sim_path, sol_path; t = t, n_itp = config.n_itp)

    summarize_result(stdout, result, connected; label = "TimCop 64^3 validation plot")
    println("  time index: $t")
    output_path = save_detection_figure(
        psi,
        X,
        result,
        config.output_path;
        isovalue = 0.5,
        isorange = 0.12,
    )
    println("  wrote plot: $output_path")
end

main(ARGS)
