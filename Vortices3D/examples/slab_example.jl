using VortexDistributions
using JLD2
using GLMakie

@load "data/slab_tubes.jld2" psi_tubes1 psi_tubes2 psi_tubes3 X 

include("../util.jl")

psi = psi_tubes1;
# psi = psi_tubes2;
# psi = psi_tubes3;

plot_iso(psi)

vortsx = findvortices3D_x(psi, X);
vortsy = findvortices3D_y(psi, X);
vortsz = findvortices3D_z(psi, X);

plot_vfound3D_box(vortsx, X, true)
plot_vfound3D_box(vortsy, X, true)
plot_vfound3D_box(vortsz, X, true)