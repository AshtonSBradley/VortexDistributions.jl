using Interpolations, FLoops, Graphs

function vortex_coords(vorts_x, vorts_y, vorts_z, x, y, z)
    dx = x[2] - x[1]; dy = y[2] - y[1]; dz = z[2] - z[1];
    vorts_x_coords = [[x[v[1]], y[v[2]] - dy/2 , z[v[3]] - dz/2] for v in vorts_x]
    vorts_y_coords = [[x[v[1]] - dx/2, y[v[2]], z[v[3]] - dz/2] for v in vorts_y]
    vorts_z_coords = [[x[v[1]] - dx/2, y[v[2]] - dy/2, z[v[3]]] for v in vorts_z]
    vorts_coords = vcat(vorts_x_coords, vorts_y_coords, vorts_z_coords);
    return vorts_coords
end

function link_graph_vorts(g; repeat_end_of_ring=true)
    g_temp = deepcopy(g)
    visited = Set()
    deg_1 = Set(findall(x -> x == 1, degree(g_temp)))

    # Link vortex lines by starting at degree 1 vortices (end points) and terminating at degree != 2 vertices
    vort_lines = []
    while length(deg_1) > 0
        current_vort = []
        vc = pop!(deg_1)
        
        push!(current_vort, vc)
        vc_neighbors = neighbors(g_temp, vc)

        # setdiff!(vc_neighbors, visited)
        while length(vc_neighbors) <= 2 # vc isn't a reconnection
            push!(visited, vc)
            setdiff!(vc_neighbors, visited)
            # print(length(vc_neighbors))
            if length(vc_neighbors) == 0
                break
            end
            vc = pop!(vc_neighbors)
            push!(current_vort, vc)
            vc_neighbors = neighbors(g_temp, vc)  
        end
        push!(vort_lines, current_vort)
        setdiff!(deg_1, visited)
    end
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

function vortex_edge_list_threaded(Δϕx, Δϕy, Δϕz)
    vorts_x = Tuple.(findall(Δϕx .> 0.0)) .|> collect;
    vorts_y = Tuple.(findall(Δϕy .> 0.0)) .|> collect;
    vorts_z = Tuple.(findall(Δϕz .> 0.0)) .|> collect;

    vorts_x_len = length(vorts_x)
    vorts_y_len = length(vorts_y)
    vorts_z_len = length(vorts_z)

    vorts_x_map = Dict((vorts_x[i], i) for i in eachindex(vorts_x))
    vorts_y_map = Dict((vorts_y[i], i + vorts_x_len) for i in eachindex(vorts_y))
    vorts_z_map = Dict((vorts_z[i], i + vorts_x_len + vorts_y_len) for i in eachindex(vorts_z))

    edge_list_x = [[] for _ in 1:Threads.nthreads()]
    edge_list_periodic_x = [[] for _ in 1:Threads.nthreads()]
 
    @floop for i in eachindex(vorts_x)
        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_x, i, Δϕx, vorts_x_map, [-1, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_x[Threads.threadid()], (i, v_idx)) : push!(edge_list_x[Threads.threadid()], (i, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_x, i, Δϕy, vorts_y_map, [0, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_x[Threads.threadid()], (i, v_idx)) : push!(edge_list_x[Threads.threadid()], (i, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_x, i, Δϕy, vorts_y_map, [0, -1, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_x[Threads.threadid()], (i, v_idx)) : push!(edge_list_x[Threads.threadid()], (i, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_x, i, Δϕz, vorts_z_map, [0, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_x[Threads.threadid()], (i, v_idx)) : push!(edge_list_x[Threads.threadid()], (i, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_x, i, Δϕz, vorts_z_map, [0, 0, -1])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_x[Threads.threadid()], (i, v_idx)) : push!(edge_list_x[Threads.threadid()], (i, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_x, i, Δϕx, vorts_x_map, [1, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_x[Threads.threadid()], (i, v_idx)) : push!(edge_list_x[Threads.threadid()], (i, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_x, i, Δϕy, vorts_y_map, [1, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_x[Threads.threadid()], (i, v_idx)) : push!(edge_list_x[Threads.threadid()], (i, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_x, i, Δϕy, vorts_y_map, [1, -1, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_x[Threads.threadid()], (i, v_idx)) : push!(edge_list_x[Threads.threadid()], (i, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_x, i, Δϕz, vorts_z_map, [1, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_x[Threads.threadid()], (i, v_idx)) : push!(edge_list_x[Threads.threadid()], (i, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_x, i, Δϕz, vorts_z_map, [1, 0, -1])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_x[Threads.threadid()], (i, v_idx)) : push!(edge_list_x[Threads.threadid()], (i, v_idx)))
    end

    edge_list_y = [[] for _ in 1:Threads.nthreads()]
    edge_list_periodic_y = [[] for _ in 1:Threads.nthreads()]

    @floop for i in eachindex(vorts_y)
        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_y, i, Δϕy, vorts_y_map, [0, -1, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_y[Threads.threadid()], (i + vorts_x_len, v_idx)) : push!(edge_list_y[Threads.threadid()], (i + vorts_x_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_y, i, Δϕz, vorts_z_map, [0, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_y[Threads.threadid()], (i + vorts_x_len, v_idx)) : push!(edge_list_y[Threads.threadid()], (i + vorts_x_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_y, i, Δϕz, vorts_z_map, [0, 0, -1])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_y[Threads.threadid()], (i + vorts_x_len, v_idx)) : push!(edge_list_y[Threads.threadid()], (i + vorts_x_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_y, i, Δϕx, vorts_x_map, [0, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_y[Threads.threadid()], (i + vorts_x_len, v_idx)) : push!(edge_list_y[Threads.threadid()], (i + vorts_x_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_y, i, Δϕx, vorts_x_map, [-1, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_y[Threads.threadid()], (i + vorts_x_len, v_idx)) : push!(edge_list_y[Threads.threadid()], (i + vorts_x_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_y, i, Δϕy, vorts_y_map, [0, 1, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_y[Threads.threadid()], (i + vorts_x_len, v_idx)) : push!(edge_list_y[Threads.threadid()], (i + vorts_x_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_y, i, Δϕz, vorts_z_map, [0, 1, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_y[Threads.threadid()], (i + vorts_x_len, v_idx)) : push!(edge_list_y[Threads.threadid()], (i + vorts_x_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_y, i, Δϕz, vorts_z_map, [0, 1, -1])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_y[Threads.threadid()], (i + vorts_x_len, v_idx)) : push!(edge_list_y[Threads.threadid()], (i + vorts_x_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_y, i, Δϕx, vorts_x_map, [0, 1, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_y[Threads.threadid()], (i + vorts_x_len, v_idx)) : push!(edge_list_y[Threads.threadid()], (i + vorts_x_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_y, i, Δϕx, vorts_x_map, [-1, 1, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_y[Threads.threadid()], (i + vorts_x_len, v_idx)) : push!(edge_list_y[Threads.threadid()], (i + vorts_x_len, v_idx)))

    end

    edge_list_z = [[] for _ in 1:Threads.nthreads()]
    edge_list_periodic_z = [[] for _ in 1:Threads.nthreads()]

    @floop for i in eachindex(vorts_z)
        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_z, i, Δϕz, vorts_z_map, [0, 0, -1])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)) : push!(edge_list_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_z, i, Δϕx, vorts_x_map, [0, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)) : push!(edge_list_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_z, i, Δϕx, vorts_x_map, [-1, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)) : push!(edge_list_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_z, i, Δϕy, vorts_y_map, [0, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)) : push!(edge_list_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_z, i, Δϕy, vorts_y_map, [0, -1, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)) : push!(edge_list_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_z, i, Δϕz, vorts_z_map, [0, 0, 1])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)) : push!(edge_list_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_z, i, Δϕx, vorts_x_map, [0, 0, 1])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)) : push!(edge_list_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_z, i, Δϕx, vorts_x_map, [-1, 0, 1])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)) : push!(edge_list_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_z, i, Δϕy, vorts_y_map, [0, 0, 1])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)) : push!(edge_list_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_z, i, Δϕy, vorts_y_map, [0, -1, 1])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)) : push!(edge_list_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)))
    end
    
    edge_list_x = reduce(vcat, edge_list_x)
    edge_list_y = reduce(vcat, edge_list_y)
    edge_list_z = reduce(vcat, edge_list_z)
    edge_list_periodic_x = reduce(vcat, edge_list_periodic_x)
    edge_list_periodic_y = reduce(vcat, edge_list_periodic_y)
    edge_list_periodic_z = reduce(vcat, edge_list_periodic_z)

    return vcat(edge_list_x, edge_list_y, edge_list_z), vcat(edge_list_periodic_x, edge_list_periodic_y, edge_list_periodic_z), vorts_x, vorts_y, vorts_z
end

function findvortices_jumps_plane(phase)
    # phase = angle.(ψ);

    Δϕx, Δϕy = phase_jumps(phase,1),phase_jumps(phase,2)

    circshift!(phase,Δϕx,(0,1))
    Δϕx .-= phase; Δϕx .-= Δϕy
    circshift!(phase,Δϕy,(1,0))
    Δϕx .+= phase

    return abs.(Δϕx)
end

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

function full_algorithm(psi, x, y, z; n_itp = 1)
    if n_itp > 1
        range_x = collect(range(start=0, stop=length(x), length=length(x)*n_itp + 1))[2:end]
        range_y = collect(range(start=0, stop=length(y), length=length(y)*n_itp + 1))[2:end]
        range_z = collect(range(start=0, stop=length(z), length=length(z)*n_itp + 1))[2:end]

        x_itp = interpolate(x, BSpline(Linear()))
        y_itp = interpolate(y, BSpline(Linear()))
        z_itp = interpolate(z, BSpline(Linear()))

        psi_itp = interpolate(psi, BSpline(Linear()))

        x_etp = extrapolate(x_itp, Line())
        y_etp = extrapolate(y_itp, Line())
        z_etp = extrapolate(z_itp, Line())

        psi_etp = extrapolate(psi_itp, Line())

        x = x_etp(range_x); y = y_etp(range_y); z = z_etp(range_z)
        psi = psi_etp(range_x, range_y, range_z)
    end


    # @time Δϕx, Δϕy, Δϕz = find_vortices_circshift(psi, x, y, z);
    println("==============================")
    print("Finding vortex points on planes:")
    @time Δϕx, Δϕy, Δϕz = findvortices_planes_threaded(psi);

    print("Creating edge list:")
    @time edge_list, edge_list_periodic, vorts_x, vorts_y, vorts_z = vortex_edge_list_threaded(Δϕx, Δϕy, Δϕz);
    
    print("Creating graph:")
    @time g = SimpleGraph(Edge.(edge_list))

    print("Linking vortices:")
    @time vort_lines, vort_loops, vort_rings = link_graph_vorts(g, repeat_end_of_ring=true);

    print("Vort coords:")
    @time vorts_coords = vortex_coords(vorts_x, vorts_y, vorts_z, x, y, z)
    println("")
    println("Number of vortex points: $(length(vorts_coords))")
    println("Number of vortex lines: $(length(vort_lines))")
    println("Number of vortex loops: $(length(vort_loops))")
    println("Number of vortex rings: $(length(vort_rings))")
    println("")
    println("Number of vortices: $(length(vort_lines) + length(vort_loops) + length(vort_rings))")
    # @time connected_vorts = connect_vortex_ends(vorts_coords, vort_lines, vort_loops, vort_rings, X);
    println("==============================")
    return g, vort_lines, vort_loops, vort_rings, vorts_coords
end

