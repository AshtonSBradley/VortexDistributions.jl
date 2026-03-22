using VortexDistributions

include(joinpath(@__DIR__, "timcop_plot_helpers.jl"))
using .TimCopPlotHelpers

"""
Run a self-contained 3D detection demo using the bundled `test/3d/box_vorts.jld2`
fixture and save a GLMakie/GraphMakie plot.

Usage:

    julia --project=. examples/timcop_detection_demo.jl [output_png] [n_itp]

Defaults:
- `output_png = examples/timcop_detection_demo.png`
- `n_itp = 1`

This script is intentionally outside the package test suite. It demonstrates the
vendored Tim Cop 3D detector on a lightweight fixture using `GLMakie` and
`GraphMakie`, matching the visualization style of the original repository.
"""

function parse_args(args)
    output_path = joinpath(@__DIR__, "timcop_detection_demo.png")
    n_itp = 1

    positional = String[]
    for arg in args
        push!(positional, arg)
    end

    if !isempty(positional)
        output_path = positional[1]
    end
    if length(positional) >= 2
        n_itp = parse(Int, positional[2])
    end

    return DemoConfig(; output_path, n_itp)
end

function main(args)
    config = parse_args(args)
    psi, X, result, connected = detect_fixture_demo(; n_itp = config.n_itp)

    summarize_result(stdout, result, connected; label = "TimCop 3D detection demo")
    output_path = save_detection_figure(
        psi,
        X,
        result,
        config.output_path;
        isovalue = config.isovalue,
        isorange = config.isorange,
    )
    println("  wrote plot: $output_path")
    println("  fixture: $(default_fixture_path())")
end

main(ARGS)
