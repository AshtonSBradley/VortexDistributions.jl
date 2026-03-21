function recursive_cluster_algorithm(vortices::AbstractVector{PointVortex}, Lx::Real, Ly::Real)
    remaining = collect(vortices)
    dipoles = Dipole[]

    while true
        dipole_result = get_dipoles(remaining, Lx, Ly)
        isempty(dipole_result.pairs) && break
        append!(dipoles, dipole_result.dipoles)
        remaining = remaining[.!dipole_result.mask]
        length(remaining) < 2 && break
    end

    if length(remaining) < 2
        return RCAResult(dipoles, Cluster[], Cluster[], remaining)
    end

    seed_result = seed_clusters(remaining, Lx, Ly)
    positive_clusters = grow_plus_clusters(remaining, seed_result.pairs, Lx, Ly)
    negative_clusters = grow_minus_clusters(remaining, seed_result.pairs, Lx, Ly)

    clustered = falses(length(remaining))
    for cluster in (positive_clusters..., negative_clusters...)
        for (index, vortex) in pairs(remaining)
            if any(vortex === member for member in cluster.vortices)
                clustered[index] = true
            end
        end
    end

    return RCAResult(dipoles, positive_clusters, negative_clusters, remaining[.!clustered])
end

function recursive_cluster_algorithm(
    xi::AbstractVector,
    yi::AbstractVector,
    charges::AbstractVector{<:Integer},
    Lx::Real,
    Ly::Real,
)
    @assert length(xi) == length(yi) == length(charges)
    vortices = PointVortex.(Float64.(xi), Float64.(yi), Int.(charges))
    return recursive_cluster_algorithm(vortices, Lx, Ly)
end
