using VortexDistributions
using Base.Test

# write your own tests here
include("testallcharges.jl")
include("testallpositions.jl")
include("makepsi.jl")
include("checkpositions.jl")
tic()
@testset "Vortex location and charge tests " begin include("vortextests.jl") end
toc()
