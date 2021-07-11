module VortexDistributions

using Interpolations 
using JLD2
using Parameters
using SpecialFunctions
using LinearAlgebra
using ToeplitzMatrices
using SparseArrays
using FFTW
using FileIO
using ProgressMeter
using LightGraphs
using SimpleWeightedGraphs
import Plots:stroke,scatter!,plot,plot!

const Λ = 0.8249

# topology
export Field, Torus, Sphere,

# vortex groups
VortexGroup, Dipole, Cluster,

# cores
Vortex, CoreShape, Ansatz, Exact, ScalarVortex, PointVortex,

# detection 
findvortices, found_near, phase_jumps, phase_jumps!, unwrap, unwrap!, Δ,
find_where, findvortices_jumps, findvortices_grid,
remove_vortices_edge, zoom_interp, zoom_grid,
circ_mask, keep_vortices,

# construction
scalar_ansatz, vortex_array, uniform,
vortex!, dipole_phase, periodic_dipole!,
rand_charge, rand_pointvortex, rand_scalarvortex, rand_vortex, rand_vortexfield,
thetad, 

# convenient access
charge, xpos, ypos, pos,

# RCA
distances, periodic_distances, sparse_distances, 
spanning_tree, vortex_marker, plot_vortices!,
plot_cluster!, vortex_marker, grow_plus_clusters, 
grow_minus_clusers, get_dipoles, seed_clusters

include("types.jl")
include("pointvortex.jl")
include("detection.jl")
include("creation.jl")

# RCA
include("get_dipoles.jl")
include("seed_clusters.jl")
include("grow_minus_clusters.jl")
include("grow_plus_clusters.jl")
# include("get_spanning_trees.jl")
# include("get_negative_spanning_trees.jl")
# include("get_positive_spanning_trees.jl")
# include("find_smallest_spanning_trees.jl")
# include("sort_RCA_structs.jl")
# include("RCA.jl")

# utils
include("utils.jl")

@load joinpath(@__DIR__,"cores.jld2") ψi ψa 

end
