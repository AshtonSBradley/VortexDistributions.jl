__precompile__()

module VortexDistributions

#using Reexport
#@reexport using DifferentialEquations

include("findvortices.jl")
include("unwrap.jl")
include("vortexcore.jl")
include("makevortex.jl")
include("circmask.jl")
include("edgemask.jl")
include("findvortmask.jl")
include("countphasejumps.jl")
include("remove_edgevortices.jl")


export findvortices, unwrap, unwrap!, makevortex, makevortex!, vortexcore, circmask, edgemask, findvortmask, remove_edgevortices

end # module
