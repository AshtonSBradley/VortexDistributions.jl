module VortexDistributions

using Reexport
@reexport using Interpolations
@reexport using SpecialFunctions
@reexport using LinearAlgebra
@reexport using ToeplitzMatrices
@reexport using SparseArrays
@reexport using FFTW
@reexport using FileIO
@reexport using JLD2

Interpolations.interpolate(A::Array{Complex{Float64},2})=interpolate(real(A))+im*interpolate(imag(A))

const Λ = 0.8249

include("findvortices.jl")
include("detection.jl")
include("makevortex.jl")
include("findvortmask.jl")
include("remove_edgevortices.jl")
include("vortexcore.jl")
include("corezoom.jl")
include("helpers.jl")

export findvortices_grid, findvortices_interp, findvortices,
unwrap, unwrap!, phasejumps, makevortex, makevortex!,
makeallvortices!, vortexcore, gpecore, circmask, edgemask,
findvortmask, remove_edgevortices, linspace, findnonzero,
randomvortices, isinterior, checkvortexlocations, corezoom,
core_chargen, make_fastcore

end # module
