using VortexDistributions
using JLD2
using GLMakie

@load "data/box_rings.jld2" psi_ring1 psi_ring2 psi_ring3 X # make sure julia repl is in base directory 

## Plot iso surface for select ring 
include("../util.jl")

psi = psi_ring2;
plot_iso(psi)

vortsx = findvortices3D_x(psi, X);
vortsy = findvortices3D_y(psi, X);
vortsz = findvortices3D_z(psi, X);

# plot_vfound3D(vortsx, X, true)
# plot_vfound3D(vortsy, X, true)
plot_vfound3D(vortsz, X, true)