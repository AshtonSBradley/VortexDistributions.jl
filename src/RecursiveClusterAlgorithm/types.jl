struct RCAResult
    dipoles::Vector{Dipole}
    positive_clusters::Vector{Cluster}
    negative_clusters::Vector{Cluster}
    free_vortices::Vector{PointVortex}
end
