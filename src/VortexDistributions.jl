module VortexDistributions

using Reexport
@reexport using Test
@reexport using JLD2
@reexport using Parameters
@reexport using Interpolations
Core.eval(Main, :(using Interpolations))
@reexport using SpecialFunctions
@reexport using LinearAlgebra
@reexport using ToeplitzMatrices
@reexport using SparseArrays
@reexport using FFTW
@reexport using FileIO

const Λ = 0.8249
export Field, Torus, Sphere
export Vortex, CoreShape, Ansatz, Exact, ScalarVortex
export PointVortex, RawData, uniform, randcharge
export randPointVortex, randScalarVortex, randVortex
export vortex!, findvortices
export foundNear, randVortexField

# export ψi, ψa
# export findwhere, findvortices_jumps, findvortices_grid, findvortices_interp
# export removeedgevortices, corezoom
# export phasejumps, phasejumps!, unwrap, unwrap!
# export gpecore_exact, chebdif, getChebDMatrix, getChebD2Matrix

Interpolations.interpolate(A::Array{Complex{Float64},2}) = interpolate(real(A))+im*interpolate(imag(A))

include("types.jl")
include("pointvortex.jl")
include("detection.jl")
include("creation.jl")

@load joinpath(@__DIR__,"exactcore.jld2") ψi
@load joinpath(@__DIR__,"ansatzcore.jld2") ψa

end # module
