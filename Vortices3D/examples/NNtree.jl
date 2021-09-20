using VortexDistributions, JLD2, Interpolations;
using NearestNeighbors, LinearAlgebra, Distances;
using GLMakie

include("../util.jl");
@load "Vortices3D/data/box_vorts.jld2" psi_chaos1 psi_chaos2 psi_knots1 psi_knots2 psi_knots3 psi_knots4 psi_tubes1 psi_tubes2 psi_ringtube X
# @load "data/box_rings.jld2" psi_ring1 psi_ring2 psi_ring3 X # make sure julia repl is in base directory 


psi = psi_knots4;
plot_iso(psi,X)

x = X[1]; y = X[2]; z = X[3];
Δx = x[2]-x[1]; Δy = y[2]-y[1]; Δz = z[2]-z[1];

N = 2
α = 0.05
ϵ = (1+α)*sqrt(Δx^2+Δy^2+Δz^2)/N;

@time nearest = uniqueVortices(ϵ, psi, X, N)
@time sets1 = setMethod(psi, X, N, ϵ)
@time sets2 = setMethod2(psi, X, N, ϵ)
@time sets3 = setMethod3(psi, X, N, ϵ)


@time vorts_3d = findvortices3D_itp(psi, X, N) # Find vortex points with interpolation depth N
v_matrix = vcat(vorts_3d'...)[:,1:3]' # Convert to matrix for kdtree 
sortslices(v_matrix, dims=2, by=x->x[1])


