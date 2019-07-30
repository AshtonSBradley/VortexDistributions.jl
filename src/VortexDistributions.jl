module VortexDistributions

using Reexport
@reexport using Test
@reexport using Parameters
@reexport using Interpolations
@reexport using SpecialFunctions
@reexport using LinearAlgebra
@reexport using ToeplitzMatrices
@reexport using SparseArrays
@reexport using FFTW
@reexport using FileIO
@reexport using JLD2

Interpolations.interpolate(A::Array{Complex{Float64},2}) = interpolate(real(A))+im*interpolate(imag(A))
@load "./src/exactcore.jld2" ψi
@load "./src/ansatzcore.jld2" ψa

include("types.jl")
include("pointvortex.jl")
include("detection.jl")
include("creation.jl")

export Field, Torus, Sphere
export Vortex, CoreShape, Anstaz, Exact, ScalarVortex, PointVortex
export randPointVortex, randScalarVortex, randVortex, vortex!
export phasejumps, phasejumps!, unwrap, unwrap!
export PointVortex, RawData, uniform, randcharge, randPointVortex
export findwhere, findvortices_jumps, findvortices_grid, findvortices_interp
export findvortices, removeedgevortices, corezoom
export ψi, ψa, gpecore_exact, chebdif, getChebDMatrix, getChebD2Matrix

# include("findvortices.jl")
# include("detection.jl")
# include("makevortex.jl")
# include("findvortmask.jl")
# include("remove_edgevortices.jl")
# include("vortexcore.jl")
# include("corezoom.jl")
# include("helpers.jl")

# export findvortices_grid, findvortices_interp, findvortices,
# unwrap, unwrap!, phasejumps, makevortex, makevortex!,
# makeallvortices!, vortexcore, gpecore, circmask, edgemask,
# findvortmask, remove_edgevortices, linspace, findnonzero,
# randomvortices, isinterior, checkvortexlocations, corezoom,
# core_chargen, make_fastcore

end # module
