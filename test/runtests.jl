using VortexDistributions
using Base.Test

# write your own tests here
include("testall_charges.jl")
include("testall_positions.jl")
include("makepsi.jl")
include("checkvortexlocations.jl")
@testset "Vortex location and circulation tests " begin include("vortextests.jl") end
