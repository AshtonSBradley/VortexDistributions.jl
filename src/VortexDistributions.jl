module VortexDistributions

using Test
using JLD2
using Parameters
using SpecialFunctions
using Interpolations
using LinearAlgebra
using ToeplitzMatrices
using SparseArrays
using FFTW
using FileIO
using ProgressMeter

const Λ = 0.8249
export Field, Torus, Sphere
export Vortex, CoreShape, Ansatz, Exact, ScalarVortex, PointVortex
export scalar_ansatz, vortex_array, uniform
export vortex!, findvortices, dipole_phase, periodic_dipole!
export rand_charge, rand_pointvortex, rand_scalarvortex, rand_vortex, rand_vortexfield
export found_near, phase_jumps, phase_jumps!, unwrap, unwrap!, Δ
export Dipole, Cluster, charge, xpos, ypos, pos

# export ψi, ψa
# export findwhere, findvortices_jumps, findvortices_grid, findvortices_interp
# export removeedgevortices, corezoom
# export gpecore_exact, chebdif, getChebDMatrix, getChebD2Matrix

include("types.jl")
include("pointvortex.jl")
include("detection.jl")
include("creation.jl")
include("utils.jl")

@load joinpath(@__DIR__,"exactcore.jld2") ψi
@load joinpath(@__DIR__,"ansatzcore.jld2") ψa

end
