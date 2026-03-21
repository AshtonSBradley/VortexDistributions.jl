# Vendored from timcop/VortexDetection.jl (commit 5ce3982), trimmed to the
# package-core 3D detection path and adapted to the VortexDistributions layout.

function findvortices_planes_threaded(ψ; n_itr = 1)
    ϕ = angle.(ψ)

    Δϕx = zeros(size(ψ))
    Δϕy = zeros(size(ψ))
    Δϕz = zeros(size(ψ))

    @floop for i in eachindex(ψ[:, 1, 1])
        Δϕx[i, :, :] = findvortices_jumps_plane(ϕ[i, :, :])
    end

    @floop for j in eachindex(ψ[1, :, 1])
        Δϕy[:, j, :] = findvortices_jumps_plane(ϕ[:, j, :])
    end

    @floop for k in eachindex(ψ[1, 1, :])
        Δϕz[:, :, k] = findvortices_jumps_plane(ϕ[:, :, k])
    end

    return Δϕx, Δϕy, Δϕz
end

function link_graph_vorts(g; repeat_end_of_ring = true)
    g_temp = deepcopy(g)
    visited = Set()
    deg_1 = Set(findall(x -> x == 1, degree(g_temp)))

    vort_lines = []
    while length(deg_1) > 0
        current_vort = []
        vc = pop!(deg_1)

        push!(current_vort, vc)
        vc_neighbors = Set(neighbors(g_temp, vc))

        while length(vc_neighbors) <= 2
            push!(visited, vc)
            setdiff!(vc_neighbors, visited)
            if length(vc_neighbors) == 0
                break
            end
            vc = pop!(vc_neighbors)
            push!(current_vort, vc)
            vc_neighbors = Set(neighbors(g_temp, vc))
        end
        push!(vort_lines, current_vort)
        setdiff!(deg_1, visited)
    end

    g_temp = deepcopy(g)
    deg_2g = Set(findall(x -> x > 2, degree(g_temp)))
    vort_loops = []

    while length(deg_2g) > 0
        current_vort = []
        vi = pop!(deg_2g)
        push!(current_vort, vi)
        push!(visited, vi)

        vi_neighbors = Set(neighbors(g_temp, vi))

        vi = nothing
        for v in vi_neighbors
            if v ∉ visited && length(neighbors(g_temp, v)) <= 2
                vi = v
                break
            end
        end

        while !isnothing(vi)
            push!(current_vort, vi)
            push!(visited, vi)
            vi_neighbors = Set(neighbors(g_temp, vi))
            setdiff!(vi_neighbors, visited)
            if length(vi_neighbors) == 1
                vi = pop!(vi_neighbors)
                if length(neighbors(g_temp, vi)) > 2
                    push!(current_vort, vi)
                    push!(visited, vi)
                    vi = nothing
                end
            else
                vi = nothing
            end
        end

        push!(vort_loops, current_vort)
        setdiff!(deg_2g, visited)
    end

    g_temp = deepcopy(g)

    deg_2 = Set(findall(x -> x == 2, degree(g_temp)))
    setdiff!(deg_2, visited)
    vort_rings = []
    while length(deg_2) > 0
        current_vort = []
        vi = pop!(deg_2)
        push!(current_vort, vi)
        push!(visited, vi)
        vi_neighbors = Set(neighbors(g_temp, vi))
        setdiff!(vi_neighbors, visited)
        while length(vi_neighbors) != 0
            vc = pop!(vi_neighbors)
            push!(current_vort, vc)
            push!(visited, vc)
            vi_neighbors = Set(neighbors(g_temp, vc))
            setdiff!(vi_neighbors, visited)
        end
        push!(vort_rings, current_vort)
        if repeat_end_of_ring
            push!(current_vort, current_vort[1])
        end
        deg_2 = setdiff!(deg_2, visited)
    end

    return vort_lines, vort_loops, vort_rings
end

function connect_vortex_ends(vorts, vort_lines, vort_loops, vort_rings, X)
    x, y, z = X
    dx = x[2] - x[1]
    dy = y[2] - y[1]
    dz = z[2] - z[1]
    Lx = (x[end] - x[1]) + dx
    Ly = (y[end] - y[1]) + dy
    Lz = (z[end] - z[1]) + dz

    vorts_lines_loops = vcat(vort_lines, vort_loops)
    line_ends = [[vorts[v[1]], vorts[v[end]]] for v in vorts_lines_loops]

    connection_array = zeros(Bool, length(vorts_lines_loops), length(vorts_lines_loops))
    for i in eachindex(vorts_lines_loops)
        for j in 1:2
            vij = line_ends[i][j]

            is_x_vort = vij[1] % dx ≈ 0
            is_y_vort = vij[2] % dy ≈ 0
            is_z_vort = vij[3] % dz ≈ 0

            x_dist_1 = is_x_vort ? dx : dx / 2
            y_dist_1 = is_y_vort ? dy : dy / 2
            z_dist_1 = is_z_vort ? dz : dz / 2

            for ii in eachindex(vorts_lines_loops)
                for jj in 1:2
                    vjj = line_ends[ii][jj]

                    x_dist = abs(vij[1] - vjj[1])
                    y_dist = abs(vij[2] - vjj[2])
                    z_dist = abs(vij[3] - vjj[3])

                    if (x_dist <= x_dist_1 || x_dist >= Lx - x_dist_1) &&
                       (y_dist <= y_dist_1 || y_dist >= Ly - y_dist_1) &&
                       (z_dist <= z_dist_1 || z_dist >= Lz - z_dist_1)
                        connection_array[i, ii] = true
                    end
                end
            end
        end
    end

    for i in axes(connection_array, 1)
        connection_array[i, i] = false
    end

    connection_graph = SimpleGraph(connection_array)
    connected_vorts_idxs = connected_components(connection_graph)

    connected_vorts = []
    for idxs in connected_vorts_idxs
        push!(connected_vorts, vorts_lines_loops[idxs])
    end

    for ring in vort_rings
        push!(connected_vorts, [ring])
    end

    return connected_vorts
end

function findvortices_jumps_plane(phase)
    Δϕx, Δϕy = phase_jumps(phase, 1), phase_jumps(phase, 2)

    circshift!(phase, Δϕx, (0, 1))
    Δϕx .-= phase
    Δϕx .-= Δϕy
    circshift!(phase, Δϕy, (1, 0))
    Δϕx .+= phase

    return abs.(Δϕx)
end

function vortex_coords(vorts_x, vorts_y, vorts_z, x, y, z)
    dx = x[2] - x[1]
    dy = y[2] - y[1]
    dz = z[2] - z[1]
    vorts_x_coords = [[x[v[1]], y[v[2]] - dy / 2, z[v[3]] - dz / 2] for v in vorts_x]
    vorts_y_coords = [[x[v[1]] - dx / 2, y[v[2]], z[v[3]] - dz / 2] for v in vorts_y]
    vorts_z_coords = [[x[v[1]] - dx / 2, y[v[2]] - dy / 2, z[v[3]]] for v in vorts_z]
    return vcat(vorts_x_coords, vorts_y_coords, vorts_z_coords)
end

function neighbour_vort(vorts, i, Δϕd, vorts_map_d, Δneighbour)
    v_idx = vorts[i]
    v_idx_i = v_idx .+ Δneighbour
    v_idx_i_mod = mod1.(v_idx_i, size(Δϕd))
    if Δϕd[v_idx_i_mod[1], v_idx_i_mod[2], v_idx_i_mod[3]] > 0.0
        if v_idx_i == v_idx_i_mod
            return true, false, vorts_map_d[v_idx_i_mod]
        else
            return true, true, vorts_map_d[v_idx_i_mod]
        end
    else
        return false, false, nothing
    end
end

function vortex_edge_list_threaded(Δϕx::AbstractArray{<:Real, 3},
                                   Δϕy::AbstractArray{<:Real, 3},
                                   Δϕz::AbstractArray{<:Real, 3})
    vorts_x = collect(Tuple.(findall(Δϕx .> 0.0)))
    vorts_y = collect(Tuple.(findall(Δϕy .> 0.0)))
    vorts_z = collect(Tuple.(findall(Δϕz .> 0.0)))

    vx = length(vorts_x)
    vy = length(vorts_y)

    vorts_x_map = Dict((vorts_x[i], i) for i in eachindex(vorts_x))
    vorts_y_map = Dict((vorts_y[i], i + vx) for i in eachindex(vorts_y))
    vorts_z_map = Dict((vorts_z[i], i + vx + vy) for i in eachindex(vorts_z))

    nth_max = Threads.maxthreadid()

    edge_list_x = [Vector{Tuple{Int, Int}}() for _ in 1:nth_max]
    edge_list_periodic_x = [Vector{Tuple{Int, Int}}() for _ in 1:nth_max]
    edge_list_y = [Vector{Tuple{Int, Int}}() for _ in 1:nth_max]
    edge_list_periodic_y = [Vector{Tuple{Int, Int}}() for _ in 1:nth_max]
    edge_list_z = [Vector{Tuple{Int, Int}}() for _ in 1:nth_max]
    edge_list_periodic_z = [Vector{Tuple{Int, Int}}() for _ in 1:nth_max]

    @threads :static for i in eachindex(vorts_x)
        tid = threadid()
        isN, isP, v = neighbour_vort(vorts_x, i, Δϕx, vorts_x_map, (-1, 0, 0))
        isN && (isP ? push!(edge_list_periodic_x[tid], (i, v)) : push!(edge_list_x[tid], (i, v)))

        isN, isP, v = neighbour_vort(vorts_x, i, Δϕy, vorts_y_map, (0, 0, 0))
        isN && (isP ? push!(edge_list_periodic_x[tid], (i, v)) : push!(edge_list_x[tid], (i, v)))

        isN, isP, v = neighbour_vort(vorts_x, i, Δϕy, vorts_y_map, (0, -1, 0))
        isN && (isP ? push!(edge_list_periodic_x[tid], (i, v)) : push!(edge_list_x[tid], (i, v)))

        isN, isP, v = neighbour_vort(vorts_x, i, Δϕz, vorts_z_map, (0, 0, 0))
        isN && (isP ? push!(edge_list_periodic_x[tid], (i, v)) : push!(edge_list_x[tid], (i, v)))

        isN, isP, v = neighbour_vort(vorts_x, i, Δϕz, vorts_z_map, (0, 0, -1))
        isN && (isP ? push!(edge_list_periodic_x[tid], (i, v)) : push!(edge_list_x[tid], (i, v)))

        isN, isP, v = neighbour_vort(vorts_x, i, Δϕx, vorts_x_map, (1, 0, 0))
        isN && (isP ? push!(edge_list_periodic_x[tid], (i, v)) : push!(edge_list_x[tid], (i, v)))

        isN, isP, v = neighbour_vort(vorts_x, i, Δϕy, vorts_y_map, (1, 0, 0))
        isN && (isP ? push!(edge_list_periodic_x[tid], (i, v)) : push!(edge_list_x[tid], (i, v)))

        isN, isP, v = neighbour_vort(vorts_x, i, Δϕy, vorts_y_map, (1, -1, 0))
        isN && (isP ? push!(edge_list_periodic_x[tid], (i, v)) : push!(edge_list_x[tid], (i, v)))

        isN, isP, v = neighbour_vort(vorts_x, i, Δϕz, vorts_z_map, (1, 0, 0))
        isN && (isP ? push!(edge_list_periodic_x[tid], (i, v)) : push!(edge_list_x[tid], (i, v)))

        isN, isP, v = neighbour_vort(vorts_x, i, Δϕz, vorts_z_map, (1, 0, -1))
        isN && (isP ? push!(edge_list_periodic_x[tid], (i, v)) : push!(edge_list_x[tid], (i, v)))
    end

    @threads :static for i in eachindex(vorts_y)
        tid = threadid()
        base = i + vx

        isN, isP, v = neighbour_vort(vorts_y, i, Δϕy, vorts_y_map, (0, -1, 0))
        isN && (isP ? push!(edge_list_periodic_y[tid], (base, v)) : push!(edge_list_y[tid], (base, v)))

        isN, isP, v = neighbour_vort(vorts_y, i, Δϕz, vorts_z_map, (0, 0, 0))
        isN && (isP ? push!(edge_list_periodic_y[tid], (base, v)) : push!(edge_list_y[tid], (base, v)))

        isN, isP, v = neighbour_vort(vorts_y, i, Δϕz, vorts_z_map, (0, 0, -1))
        isN && (isP ? push!(edge_list_periodic_y[tid], (base, v)) : push!(edge_list_y[tid], (base, v)))

        isN, isP, v = neighbour_vort(vorts_y, i, Δϕx, vorts_x_map, (0, 0, 0))
        isN && (isP ? push!(edge_list_periodic_y[tid], (base, v)) : push!(edge_list_y[tid], (base, v)))

        isN, isP, v = neighbour_vort(vorts_y, i, Δϕx, vorts_x_map, (-1, 0, 0))
        isN && (isP ? push!(edge_list_periodic_y[tid], (base, v)) : push!(edge_list_y[tid], (base, v)))

        isN, isP, v = neighbour_vort(vorts_y, i, Δϕy, vorts_y_map, (0, 1, 0))
        isN && (isP ? push!(edge_list_periodic_y[tid], (base, v)) : push!(edge_list_y[tid], (base, v)))

        isN, isP, v = neighbour_vort(vorts_y, i, Δϕz, vorts_z_map, (0, 1, 0))
        isN && (isP ? push!(edge_list_periodic_y[tid], (base, v)) : push!(edge_list_y[tid], (base, v)))

        isN, isP, v = neighbour_vort(vorts_y, i, Δϕz, vorts_z_map, (0, 1, -1))
        isN && (isP ? push!(edge_list_periodic_y[tid], (base, v)) : push!(edge_list_y[tid], (base, v)))

        isN, isP, v = neighbour_vort(vorts_y, i, Δϕx, vorts_x_map, (0, 1, 0))
        isN && (isP ? push!(edge_list_periodic_y[tid], (base, v)) : push!(edge_list_y[tid], (base, v)))

        isN, isP, v = neighbour_vort(vorts_y, i, Δϕx, vorts_x_map, (-1, 1, 0))
        isN && (isP ? push!(edge_list_periodic_y[tid], (base, v)) : push!(edge_list_y[tid], (base, v)))
    end

    @threads :static for i in eachindex(vorts_z)
        tid = threadid()
        base = i + vx + vy

        isN, isP, v = neighbour_vort(vorts_z, i, Δϕz, vorts_z_map, (0, 0, -1))
        isN && (isP ? push!(edge_list_periodic_z[tid], (base, v)) : push!(edge_list_z[tid], (base, v)))

        isN, isP, v = neighbour_vort(vorts_z, i, Δϕx, vorts_x_map, (0, 0, 0))
        isN && (isP ? push!(edge_list_periodic_z[tid], (base, v)) : push!(edge_list_z[tid], (base, v)))

        isN, isP, v = neighbour_vort(vorts_z, i, Δϕx, vorts_x_map, (-1, 0, 0))
        isN && (isP ? push!(edge_list_periodic_z[tid], (base, v)) : push!(edge_list_z[tid], (base, v)))

        isN, isP, v = neighbour_vort(vorts_z, i, Δϕy, vorts_y_map, (0, 0, 0))
        isN && (isP ? push!(edge_list_periodic_z[tid], (base, v)) : push!(edge_list_z[tid], (base, v)))

        isN, isP, v = neighbour_vort(vorts_z, i, Δϕy, vorts_y_map, (0, -1, 0))
        isN && (isP ? push!(edge_list_periodic_z[tid], (base, v)) : push!(edge_list_z[tid], (base, v)))

        isN, isP, v = neighbour_vort(vorts_z, i, Δϕz, vorts_z_map, (0, 0, 1))
        isN && (isP ? push!(edge_list_periodic_z[tid], (base, v)) : push!(edge_list_z[tid], (base, v)))

        isN, isP, v = neighbour_vort(vorts_z, i, Δϕx, vorts_x_map, (0, 0, 1))
        isN && (isP ? push!(edge_list_periodic_z[tid], (base, v)) : push!(edge_list_z[tid], (base, v)))

        isN, isP, v = neighbour_vort(vorts_z, i, Δϕx, vorts_x_map, (-1, 0, 1))
        isN && (isP ? push!(edge_list_periodic_z[tid], (base, v)) : push!(edge_list_z[tid], (base, v)))

        isN, isP, v = neighbour_vort(vorts_z, i, Δϕy, vorts_y_map, (0, 0, 1))
        isN && (isP ? push!(edge_list_periodic_z[tid], (base, v)) : push!(edge_list_z[tid], (base, v)))

        isN, isP, v = neighbour_vort(vorts_z, i, Δϕy, vorts_y_map, (0, -1, 1))
        isN && (isP ? push!(edge_list_periodic_z[tid], (base, v)) : push!(edge_list_z[tid], (base, v)))
    end

    edge_list_x = reduce(vcat, edge_list_x; init = Tuple{Int, Int}[])
    edge_list_y = reduce(vcat, edge_list_y; init = Tuple{Int, Int}[])
    edge_list_z = reduce(vcat, edge_list_z; init = Tuple{Int, Int}[])
    edge_list_periodic_x = reduce(vcat, edge_list_periodic_x; init = Tuple{Int, Int}[])
    edge_list_periodic_y = reduce(vcat, edge_list_periodic_y; init = Tuple{Int, Int}[])
    edge_list_periodic_z = reduce(vcat, edge_list_periodic_z; init = Tuple{Int, Int}[])

    return vcat(edge_list_x, edge_list_y, edge_list_z),
        vcat(edge_list_periodic_x, edge_list_periodic_y, edge_list_periodic_z),
        vorts_x,
        vorts_y,
        vorts_z
end

function _interpolate_grid(psi, x, y, z, n_itp)
    range_x = collect(range(start = 0, stop = length(x), length = length(x) * n_itp + 1))[2:end]
    range_y = collect(range(start = 0, stop = length(y), length = length(y) * n_itp + 1))[2:end]
    range_z = collect(range(start = 0, stop = length(z), length = length(z) * n_itp + 1))[2:end]

    x_itp = interpolate(x, BSpline(Linear()))
    y_itp = interpolate(y, BSpline(Linear()))
    z_itp = interpolate(z, BSpline(Linear()))
    psi_itp = interpolate(psi, BSpline(Linear()))

    x_etp = extrapolate(x_itp, Line())
    y_etp = extrapolate(y_itp, Line())
    z_etp = extrapolate(z_itp, Line())
    psi_etp = extrapolate(psi_itp, Line())

    return x_etp(range_x), y_etp(range_y), z_etp(range_z), psi_etp(range_x, range_y, range_z)
end

function full_algorithm(psi, x, y, z; n_itp = 1, repeat_end_of_ring = true)
    if n_itp > 1
        x, y, z, psi = _interpolate_grid(psi, x, y, z, n_itp)
    end

    Δϕx, Δϕy, Δϕz = findvortices_planes_threaded(psi)
    edge_list, edge_list_periodic, vorts_x, vorts_y, vorts_z = vortex_edge_list_threaded(Δϕx, Δϕy, Δϕz)
    g = SimpleGraph(Edge.(edge_list))
    vort_lines, vort_loops, vort_rings = link_graph_vorts(g, repeat_end_of_ring = repeat_end_of_ring)
    vorts_coords = vortex_coords(vorts_x, vorts_y, vorts_z, x, y, z)

    return g, vort_lines, vort_loops, vort_rings, vorts_coords, edge_list_periodic
end

function detect_vortices_3d(psi, x, y, z; n_itp = 1, repeat_end_of_ring = true)
    return Detection3DResult(full_algorithm(psi, x, y, z; n_itp = n_itp, repeat_end_of_ring = repeat_end_of_ring)...)
end
