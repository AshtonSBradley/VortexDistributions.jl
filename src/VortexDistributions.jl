__precompile__()

module VortexDistributions

#using Reexport
#@reexport using DifferentialEquations

include("findvortices.jl")
include("unwrap.jl")
include("makevortex.jl")
include("vortexcore.jl")

export findvortices, unwrap, makevortex, vortexcore

end # module
