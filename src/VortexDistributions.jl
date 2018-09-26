__precompile__()

module VortexDistributions

#using Reexport
#@reexport using DifferentialEquations

include("findvortices.jl")
include("unwrap.jl")
include("makevortex.jl")
include("circmask.jl")
include("findvortmask.jl")
include("countphasejumps.jl")
include("remove_edgevortices.jl")
include("helpers.jl")


export findvortices, unwrap, unwrap!, countphasejumps, makevortex, makevortex!, vortexcore, circmask, edgemask, findvortmask, remove_edgevortices, linspace

end # module
