using VortexDistributions
using Test

# write your own tests here
include("testall_charges.jl")
include("testall_positions.jl")
include("makepsi.jl")
include("checkvortexlocations.jl")
@testset "Vortex location and circulation tests " begin include("vortextests.jl") end

#=
Nv = 5
x,y,psi,testvort = makepsi(Nv)
testvort
nt,np,nn,vortices = findvortices(psi,x,y)
vortices = remove_edgevortices(vortices,x,y)
chargesfound = (vortices[:,3] == testvort[:,3])
=#
