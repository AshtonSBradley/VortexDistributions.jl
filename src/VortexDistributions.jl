module VortexDistributions

using Reexport
@reexport using Interpolations
@reexport using SpecialFunctions
@reexport using LinearAlgebra
@reexport using ToeplitzMatrices
@reexport using SparseArrays
@reexport using FFTW

Interpolations.interpolate(A::Array{Complex{Float64},2})=interpolate(real(A))+im*interpolate(imag(A))

include("findvortices.jl")
include("unwrap.jl")
include("makevortex.jl")
include("findvortmask.jl")
include("countphasejumps.jl")
include("remove_edgevortices.jl")
include("vortexcore.jl")
include("corezoom.jl")
include("helpers.jl")
# include("velocities.jl")
# include("helmholtz.jl")
# include("energydecomposition.jl")


export findvortices_grid, findvortices_interp, findvortices,
unwrap, unwrap!, countphasejumps, makevortex, makevortex!,
makeallvortices!, vortexcore, gpecore, circmask, edgemask,
findvortmask, remove_edgevortices, linspace, findnonzero,
randomvortices, isinterior, checkvortexlocations, corezoom,
core_chargen, make_fastcore

# export velocities, helmholtzdecomp, energydecomp

end # module
