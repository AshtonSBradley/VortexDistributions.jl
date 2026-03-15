"""
Generate 64^3 validation data for the experimental 3D detector using FourierGPE.

Usage:

    julia --project=. examples/generate_3d_validation_data.jl [output_dir] [tf] [Nt]

Defaults:
- output_dir = examples/3d_validation_data
- tf = 80.0
- Nt = 100

This script is intentionally kept out of the package test suite because it can
take a few minutes and requires the optional `FourierGPE` dependency.
"""

using JLD2
using VortexDistributions

try
    @eval using FourierGPE
catch err
    error("This script requires the optional FourierGPE dependency.\n" *
          "Install it in your active environment with:\n" *
          "  julia --project=. -e 'using Pkg; Pkg.add(\"FourierGPE\")'\n\n" *
          "Original load error: $(sprint(showerror, err))")
end

function main(args)
    output_dir = length(args) >= 1 ? args[1] : joinpath(@__DIR__, "3d_validation_data")
    tf = length(args) >= 2 ? parse(Float64, args[2]) : 80.0
    Nt = length(args) >= 3 ? parse(Int, args[3]) : 100

    mkpath(output_dir)

    L = (16.0, 16.0, 16.0)
    N = (64, 64, 64)
    sim = Sim(L, N)
    @unpack_Sim sim

    μ = 25.0
    γ = 0.05
    t = LinRange(0.0, tf, Nt)

    x, y, z = X
    ψi = randn(N) + im * randn(N)
    ϕi = kspace(ψi, sim)

    @pack_Sim! sim = μ, γ, t, ψi, ϕi

    println("Running FourierGPE simulation for 3D validation...")
    sol = runsim(sim)

    for i in eachindex(sol)
        sol[i] = xspace(sol[i], sim)
    end

    sim_path = joinpath(output_dir, "sim_64.jld2")
    sol_path = joinpath(output_dir, "sol_64.jld2")
    @save sim_path sim
    @save sol_path sol

    println("Wrote:")
    println("  $sim_path")
    println("  $sol_path")
end

main(ARGS)
