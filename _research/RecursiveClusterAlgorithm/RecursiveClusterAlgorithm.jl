module RecursiveClusterAlgorithm

# not currently in project but can compile locally!
using VortexDistributions

using VortexDistributions:VortexGroup, Dipole, xpos, ypos

using LinearAlgebra, SparseArrays, LightGraphs, SimpleWeightedGraphs 

mutable struct Cluster <: VortexGroup
    vortices::Array{PointVortex,1}
    tree::Array{LightGraphs.SimpleGraphs.SimpleEdge{Int64},1}
end

include("get_dipoles.jl")
include("seed_clusters.jl")
include("utils.jl")
Cluster(vort::Array{PointVortex,1}) = Cluster(vort,spanning_tree(xpos(vort),ypos(vort)))

export Cluster, get_dipoles, seed_clusters, distances, periodic_distances,
sparse_distances, spanning_tree, vortex_marker, plot_vortices!,
plot_cluster!, vortex_marker

end # module
