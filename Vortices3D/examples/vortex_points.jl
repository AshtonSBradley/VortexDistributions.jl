using VortexDistributions
using JLD2
using GLMakie
using Interpolations
using NearestNeighbors
using LinearAlgebra
using Distances
## Include util functions 
include("../util.jl")


## Load psi simulation 
# @load "data/box_rings.jld2" psi_ring1 psi_ring2 psi_ring3 X # make sure julia repl is in base directory 
@load "data/box_vorts.jld2" psi_chaos1 psi_chaos2 psi_knots1 psi_knots2 psi_knots3 psi_knots4 psi_tubes1 psi_tubes2 psi_ringtube X
# @load "data/slab_tubes.jld2" psi_tubes1 psi_tubes2 psi_tubes3 X 


## Set psi 
psi = psi_tubes1;


## Plot iso surface
plot_iso(psi, X)


## Find vortices
@time vorts_3d = findvortices3D_itp(psi, X, 2);
num_vorts = length(vcat(vorts_3d...)[:,1])
print("Number of vortices found: ", length(vcat(vorts_3d...)[:,1]))


## Scatter vortices found
scatterVortsOnIso(vorts_3d)
