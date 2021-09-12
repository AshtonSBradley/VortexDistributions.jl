using VortexDistributions
using JLD2
using GLMakie

@load "data/box_rings.jld2" psi_ring1 psi_ring2 psi_ring3 X # make sure julia repl is in base directory 

## Plot iso surface for select ring 
include("util.jl")

@time vorts_3d = findvortices3D_itp(psi_ring2, X, 4);

plot_iso(psi_ring2)
plot_vfound3D(vorts_3d, X, 750)

vcat(vorts_3d[1], vorts_3d[2], vorts_3d[3])
