function grow_plus_clusters(
    vortices::AbstractVector{PointVortex},
    seed_pairs::AbstractVector{<:Tuple},
    Lx::Real,
    Ly::Real,
)
    charges = _charges(vortices)
    positive_pairs = _cluster_seed_pairs(seed_pairs, charges, 1)
    cluster_sets = [
        _grow_cluster_indices(vortices, charges, pair, 1, Lx, Ly) for pair in positive_pairs
    ]
    return _build_clusters(vortices, _merge_overlapping_clusters(cluster_sets))
end

function grow_minus_clusters(
    vortices::AbstractVector{PointVortex},
    seed_pairs::AbstractVector{<:Tuple},
    Lx::Real,
    Ly::Real,
)
    charges = _charges(vortices)
    negative_pairs = _cluster_seed_pairs(seed_pairs, charges, -1)
    cluster_sets = [
        _grow_cluster_indices(vortices, charges, pair, -1, Lx, Ly) for pair in negative_pairs
    ]
    return _build_clusters(vortices, _merge_overlapping_clusters(cluster_sets))
end
