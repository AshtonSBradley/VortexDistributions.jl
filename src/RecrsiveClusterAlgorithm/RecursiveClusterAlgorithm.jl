module RecursiveClusterAlgorithm

# using VortexDistributions # not in project but can compile locally!
using LinearAlgebra, SparseArrays, LightGraphs, SimpleWeightedGraphs, Test, Plots

export get_dipoles, seed_clusters, distances, periodic_distances,
sparse_distances, spanning_tree, vortex_marker, plot_vortices!,
plot_cluster!, vortex_marker

include("get_dipoles.jl")
include("seed_clusters.jl")
include("utils.jl")

end # module
