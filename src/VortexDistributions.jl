__precompile__()

module VortexDistributions

#using Reexport
#@reexport using DifferentialEquations

include("findvortices.jl")
include("unwrap.jl")
include("vortexcore.jl")
include("makevortex.jl")
include("makevortex!.jl")

export findvortices, unwrap, makevortex, makevortex!, vortexcore

end # module
