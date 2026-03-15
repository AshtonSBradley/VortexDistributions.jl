Δ(x) = x[2] - x[1]

function find_where(A)
    I = findall(!iszero, A)
    v = A[I]
    ix = [I[i][1] for i in eachindex(I)]
    iy = [I[i][2] for i in eachindex(I)]
    return ix, iy, v
end

function distances(x, y)
    xij = @. x - x'
    yij = @. y - y'
    return xij, yij
end

function periodic_distances!(x, y, Lx, Ly)
    make_periodic!(x, Lx)
    make_periodic!(y, Ly)
end

function make_periodic!(x, L)
    @. x[x > L / 2] -= L
    @. x[x < -L / 2] += L
end

function sparse_distances(x, y)
    xij = @. x - x'
    yij = @. y - y'
    rij = hypot.(xij, yij) |> tril |> sparse
    return rij
end

function spanning_tree(x, y)
    rij = sparse_distances(x, y)
    g = SimpleWeightedGraph(length(rij))
    for edge in eachindex(rij)
        n1, n2 = edge[1], edge[2]
        add_edge!(g, n1, n2, rij[n1, n2])
    end
    return prim_mst(g)
end
