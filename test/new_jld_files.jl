using JLD2
using Test
using VortexDistributions

@load joinpath(dirname(pathof(VortexDistributions)), "cores.jld2") ψi ψa

@test ψi(0.25) isa Real
@test ψa(0.25) isa Real
@test ψi(0.0) ≈ 0.0 atol = 1e-12
@test ψa(0.0) ≈ 0.0 atol = 1e-12
