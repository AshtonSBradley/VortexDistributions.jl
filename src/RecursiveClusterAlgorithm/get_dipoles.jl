function get_dipoles(vortices::AbstractVector{PointVortex}, Lx::Real, Ly::Real)
    xi = xpos(vortices)
    yi = ypos(vortices)
    charges = _charges(vortices)
    result = get_dipoles(xi, yi, charges, Lx, Ly)
    dipoles = Dipole[Dipole(vortices[i], vortices[j]) for (i, j) in result.pairs]
    return (; dipoles, mask = result.mask, pairs = result.pairs, nearest_neighbors = result.nearest_neighbors)
end

function get_dipoles(
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
    mask = BitVector([charges[i] * charges[nn[i]] == -1 && nn[nn[i]] == i for i in eachindex(nn)])
    pairs = _mutual_pairs(nn, mask)
    dipole_mask = falses(n)
    for (i, j) in pairs
        dipole_mask[i] = true
        dipole_mask[j] = true
    end
    return (; pairs, mask = dipole_mask, nearest_neighbors = nn)
end
