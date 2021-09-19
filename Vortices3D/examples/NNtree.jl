using VortexDistributions, JLD2, Interpolations;
using NearestNeighbors, LinearAlgebra, Distances;

include("../util.jl");
@load "data/box_vorts.jld2" psi_chaos1 psi_chaos2 psi_knots1 psi_knots2 psi_knots3 psi_knots4 psi_tubes1 psi_tubes2 psi_ringtube X
# @load "data/box_rings.jld2" psi_ring1 psi_ring2 psi_ring3 X # make sure julia repl is in base directory 


psi = psi_tubes2;

x = X[1]; y = X[2]; z = X[3];
Δx = x[2]-x[1]; Δy = y[2]-y[1]; Δz = z[2]-z[1];

N = 1
α = 0.1
ϵ = (1+α)*sqrt(Δx^2+Δy^2+Δz^2)/N;
vf = @time setMethod(psi, X, N, ϵ)

plot_iso(psi,X)
