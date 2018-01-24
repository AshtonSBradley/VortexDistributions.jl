using VortexDistributions
using Base.Test

# write your own tests here
include("testallcharges.jl")
include("testallpositions.jl")
include("makepsi.jl")
include("checkpositions.jl")
@testset "Vortex location and circulation tests " begin include("vortextests.jl") end
