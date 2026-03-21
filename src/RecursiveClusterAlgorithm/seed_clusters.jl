function seed_clusters(vortices::AbstractVector{PointVortex}, Lx::Real, Ly::Real)
    xi = xpos(vortices)
    yi = ypos(vortices)
    charges = _charges(vortices)
    return seed_clusters(xi, yi, charges, Lx, Ly)
end

function seed_clusters(
    xi::AbstractVector,
    yi::AbstractVector,
    charges::AbstractVector{<:Integer},
    Lx::Real,
    Ly::Real,
)
    n = length(xi)
    @assert length(yi) == n
    @assert length(charges) == n

    if n == 0
        return (; pairs = Tuple{Int, Int}[], mask = falses(0), nearest_neighbors = Int[])
    end

    nn = _nearest_neighbors(xi, yi, Lx, Ly)
    mask = BitVector([charges[i] * charges[nn[i]] == 1 && nn[nn[i]] == i for i in eachindex(nn)])
    pairs = _mutual_pairs(nn, mask)
    seed_mask = falses(n)
    for (i, j) in pairs
        seed_mask[i] = true
        seed_mask[j] = true
    end
    return (; pairs, mask = seed_mask, nearest_neighbors = nn)
end
