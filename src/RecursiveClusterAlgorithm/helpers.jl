function _charges(vortices::AbstractVector{PointVortex})
    return Int[charge(vortex) for vortex in vortices]
end

function _pairwise_periodic_distances(xi::AbstractVector, yi::AbstractVector, Lx::Real, Ly::Real)
    xij, yij = distances(xi, yi)
    periodic_distances!(xij, yij, Lx, Ly)
    return hypot.(xij, yij)
end

function _nearest_neighbors(xi::AbstractVector, yi::AbstractVector, Lx::Real, Ly::Real)
    n = length(xi)
    if n == 0
        return Int[]
    end
    rij = _pairwise_periodic_distances(xi, yi, Lx, Ly)
    diag = 1e6 * max(hypot(Lx, Ly), 1.0)
    rij += Diagonal(fill(diag, n))
    return Int[argmin(@view rij[i, :]) for i in 1:n]
end

function _mutual_pairs(nn::AbstractVector{<:Integer}, predicate::AbstractVector{Bool})
    pairs = Tuple{Int, Int}[]
    seen = Set{Tuple{Int, Int}}()
    for i in eachindex(nn)
        j = nn[i]
        if predicate[i] && 1 <= j <= length(nn) && nn[j] == i
            pair = i < j ? (i, j) : (j, i)
            if !(pair in seen)
                push!(pairs, pair)
                push!(seen, pair)
            end
        end
    end
    return sort(pairs)
end

function _periodic_delta(dx::Real, L::Real)
    if dx > L / 2
        return dx - L
    elseif dx < -L / 2
        return dx + L
    end
    return dx
end

function _periodic_distance(v1::PointVortex, v2::PointVortex, Lx::Real, Ly::Real)
    dx = _periodic_delta(xpos(v1) - xpos(v2), Lx)
    dy = _periodic_delta(ypos(v1) - ypos(v2), Ly)
    return hypot(dx, dy)
end

function _minimum_distance(vortex::PointVortex, candidates::AbstractVector{PointVortex}, Lx::Real, Ly::Real)
    isempty(candidates) && return Inf
    return minimum(_periodic_distance(vortex, candidate, Lx, Ly) for candidate in candidates)
end

function _dedupe_seed_pairs(seed_pairs)
    pairs = Tuple{Int, Int}[]
    seen = Set{Tuple{Int, Int}}()
    for pair in seed_pairs
        canonical = pair[1] < pair[2] ? (pair[1], pair[2]) : (pair[2], pair[1])
        if !(canonical in seen)
            push!(pairs, canonical)
            push!(seen, canonical)
        end
    end
    return sort(pairs)
end

function _cluster_seed_pairs(seed_pairs, charges, sign)
    return Tuple{Int, Int}[
        pair for pair in _dedupe_seed_pairs(seed_pairs) if charges[pair[1]] == sign && charges[pair[2]] == sign
    ]
end

function _grow_cluster_indices(
    vortices::AbstractVector{PointVortex},
    charges::AbstractVector{<:Integer},
    seed_pair::Tuple{Int, Int},
    sign::Int,
    Lx::Real,
    Ly::Real,
)
    cluster = Set([seed_pair[1], seed_pair[2]])
    same_sign_indices = findall(==(sign), charges)
    opposite_sign_indices = findall(==(-sign), charges)

    while true
        additions = Set{Int}()
        cluster_indices = sort!(collect(cluster))
        cluster_vortices = vortices[cluster_indices]
        opposite_vortices = vortices[opposite_sign_indices]

        for cluster_index in cluster_indices
            vortex = vortices[cluster_index]
            available_same_sign = [
                idx for idx in same_sign_indices if !(idx in cluster)
            ]
            isempty(available_same_sign) && continue

            distances_to_candidates = [
                _periodic_distance(vortex, vortices[idx], Lx, Ly) for idx in available_same_sign
            ]
            nearest_position = argmin(distances_to_candidates)
            candidate_index = available_same_sign[nearest_position]
            candidate = vortices[candidate_index]
            nearest_same = distances_to_candidates[nearest_position]

            nearest_opposite_to_member = _minimum_distance(vortex, opposite_vortices, Lx, Ly)
            nearest_opposite_to_candidate = _minimum_distance(candidate, opposite_vortices, Lx, Ly)

            if nearest_same < nearest_opposite_to_member && nearest_same < nearest_opposite_to_candidate
                push!(additions, candidate_index)
            end
        end

        isempty(additions) && return sort!(collect(cluster))
        union!(cluster, additions)
        length(cluster) == length(cluster_vortices) && return sort!(collect(cluster))
    end
end

function _merge_overlapping_clusters(cluster_sets::AbstractVector{Vector{Int}})
    pending = [Set(cluster) for cluster in cluster_sets if !isempty(cluster)]
    merged = Vector{Set{Int}}()
    while !isempty(pending)
        current = popfirst!(pending)
        i = 1
        while i <= length(pending)
            if !isempty(intersect(current, pending[i]))
                union!(current, pending[i])
                deleteat!(pending, i)
            else
                i += 1
            end
        end
        push!(merged, current)
    end
    return sort!.(collect.(merged))
end

function _build_clusters(vortices::AbstractVector{PointVortex}, cluster_sets::AbstractVector{Vector{Int}})
    return Cluster[Cluster(vortices[indices]) for indices in cluster_sets if length(indices) >= 2]
end
