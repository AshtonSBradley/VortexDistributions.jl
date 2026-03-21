module RecursiveClusterAlgorithm

using LinearAlgebra: Diagonal

using ..Core: Cluster,
    Dipole,
    PointVortex,
    charge,
    distances,
    periodic_distances!,
    xpos,
    ypos

include("RecursiveClusterAlgorithm/types.jl")
include("RecursiveClusterAlgorithm/helpers.jl")
include("RecursiveClusterAlgorithm/get_dipoles.jl")
include("RecursiveClusterAlgorithm/seed_clusters.jl")
include("RecursiveClusterAlgorithm/grow_clusters.jl")
include("RecursiveClusterAlgorithm/driver.jl")

export RCAResult,
    get_dipoles,
    grow_minus_clusters,
    grow_plus_clusters,
    recursive_cluster_algorithm,
    seed_clusters

end
