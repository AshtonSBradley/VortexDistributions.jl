using VortexDistributions
using JLD2
using GLMakie
using Interpolations

@load "data/box_rings.jld2" psi_ring1 psi_ring2 psi_ring3 X # make sure julia repl is in base directory 
@load "data/box_vorts.jld2" psi_chaos1 psi_chaos2 psi_knots1 psi_knots2 psi_knots3 psi_knots4 psi_tubes1 psi_tubes2 psi_ringtube X

## Plot iso surface for select ring 
include("util.jl")

psi = psi_knots4;
plot_iso(psi)

@time vorts_3d = findvortices3D_itp(psi, X, 10);

plot_vfound3D(vorts_3d, X, 750)

vcat(vorts_3d[1], vorts_3d[2], vorts_3d[3])

