
function uniqueVortices(rtol, psi, X, N)
    vorts_all = findvortices3D_itp(psi, X, N) # Find all vortex points 
    vorts = copy(vorts_all)

    filaments = [] # Array to hold unique collection's of vortices
    current = [] # The current collection 
    
    # Pop a vortex off and add it to the current collection
    vi = pop!(vorts) 
    push!(current, vi)

    count = 0 # Used to propogate down other side of filament

    while length(vorts) > 0
        dist, idx = rtol, 0 # dist and idx are used to find the nearest neighbour, initialise dist to rtol
        
        # Loop through other vortices excluding the current one
        for j in 1:length(vorts)
            vj = vorts[j]
            c_dist = sqrt((vj[1]-vi[1])^2 + (vj[2]-vi[2])^2 + (vj[3]-vi[3])^2) # euclidian distance
            
            # If we encounter a smaller distance that our current distance, update it
            if c_dist < dist
                dist, idx = c_dist, j 
            end
        end
        if dist != rtol # Found a close neighbour 
            vi = vorts[idx]
            deleteat!(vorts, idx) # Remove this vortex from vorts
            push!(current, vi) # Add it to current collection
        else  # No close neighbour found
            if count == 0 # compute the other side of vortex
                vi = current[1] # Return to starting vortex
                count+=1
            else # computed both sides, continue to next filament
                push!(filaments, current)
                current = []
                vi = pop!(vorts)
                push!(current, vi)
                count = 0
            end
        end
    end
    push!(filaments, current)
    return filaments
end


function scatterDemo(fils)
    f = copy(fils)
    color = "red"
    count = 0
    while length(f) > 0
        fi = pop!(f)
        fi_c = copy(fi)
        while length(fi_c) > 0
            vi = popfirst!(fi_c)
            scatter!([vi[1]],[vi[2]],[vi[3]],markersize=350,color=color)
            sleep(0.003)
        end
        sleep(1)
        count += 1
        if count % 2 == 0
            color = "red"
        else 
            color = "blue"
        end
    end
end

function skip_v(x)
    if x in found
        return true
    end
    return false
end

# Find all vortices within an ϵ-ball around a vortex and put this into a set.
# Then interate over the sets performing union operations if there is an intersection
function setMethod(psi, X, N, ϵ)
    vorts_3d = findvortices3D_itp(psi, X, N) # Find vortex points with interpolation depth N
    v_matrix = vcat(vorts_3d'...)[:,1:3]' # Convert to matrix for kdtree 
    num_vorts = length(vorts_3d)
    kdtree = KDTree(v_matrix)

    unvisited = Set(1:num_vorts) # Set of all unvisited vortices
    vfound = [] # Will contain sets of vortices within the ϵ-balls
    while length(unvisited) > 0
        idx = first(unvisited)
        vi = v_matrix[:, idx]
        v_idxs = inrange(kdtree, vi, ϵ) # Returns indices of vortices within radius ϵ
        if length(v_idxs) > 1 # Will always return itself, so throw away rouge vortices
            v_set = Set(v_idxs)
            push!(vfound, v_set)
        end
        setdiff!(unvisited, Set(idx)) # Remove this index from unvisited 
    end

    # A loop that unions sets that have an intersection.
    i = 1
    while i <= length(vfound)
        for j in i+1:length(vfound)
            if length(intersect(vfound[i], vfound[j])) > 0
                union!(vfound[i], vfound[j])
                deleteat!(vfound, j)
                i += -1
                break
            end
        end
        i+=1
    end
    return vfound
end

function setMethod2(v_matrix, ϵ)
    kdtree = KDTree(v_matrix)
    
    balls = inrange(kdtree, v_matrix, ϵ)
    balls = [Set(balls[i]) for i=1:length(balls) if length(balls[i]) > 1]
    # A loop that unions sets that have an intersection.
    i = 1
    while i <= length(balls)
        for j in i+1:length(balls)
            if length(intersect(balls[i], balls[j])) > 0
                union!(balls[i], balls[j])
                deleteat!(balls, j)
                i += -1
                break
            end
        end
        i+=1
    end
    return balls
end

function setMethod3(vorts_3d, ϵ)
    v_matrix = vcat(vorts_3d'...)[:,1:3]' # Convert to matrix for kdtree 
    kdtree = KDTree(v_matrix)

    # BFS using sets 
    num_vorts = length(v_matrix[1,:])
    unvisited = Set(collect(1:num_vorts))
    fils = []
    while length(unvisited) > 0
        idx = first(unvisited)
        vc = v_matrix[:, idx]
        f_idxs = inrange(kdtree, vc, ϵ)
        f = Set(f_idxs)
        search = Set(f_idxs)
        setdiff!(search, idx)
        while length(search) > 0
            idx = first(search)
            setdiff!(search, idx)
            vc = v_matrix[:, idx]
            vc_idxs = inrange(kdtree, vc, ϵ)
            setdiff!(vc_idxs, f)
            union!(f, Set(vc_idxs))
            union!(search, Set(vc_idxs))
        end
        if length(f) > 1
            push!(fils, f)
        end
        setdiff!(unvisited, f)
    end
    return fils
end


function sort_classified_vorts(v_class, vorts_3d, X)
    x = X[1]; y = X[2]; z = X[3];
    dx = x[2]-x[1]; dy = y[2]-y[1]; dz = z[2]-z[1];
    v_matrix = vcat(vorts_3d'...)[:,1:3]'

    vorts_sorted = []

    for i in 1:length(v_class)
        vi = v_class[i]
        vi = v_matrix[:, collect(vi)]
        vivecs = collect(eachcol(vi))
        visorted = []

        vc = pop!(vivecs)
        v1 = vc
        push!(visorted, vc)
        min_dist = sqrt(3)*(x[end]-x[1]) # init to max distance 
        min_idx = -1


        while length(vivecs) > 0
            min_dist = 10000
            min_idx = -1
            ##periodic
            for i in 1:length(vivecs)
                d = euclid(vc, vivecs[i]);
                if d < min_dist
                    min_dist = d
                    min_idx = i
                end
            end
            if min_dist > 2*dx
                #    if (vc[1] < x[1] || vc[1] > x[end]) || (vc[2] < y[1] || vc[2] > y[end]) || (vc[3] < z[1] || vc[3] > z[end])
                if vc[1] <= x[1]
                    vc[1] += (x[end]-x[1])
                elseif vc[1] >= x[end]
                    vc[1] += (x[1]-x[end])
                end
                if (vc[2] <= y[1])
                    vc[2] += (y[end]-y[1])

                elseif vc[2] >= y[end]
                    vc[2] += (y[1]-y[end])
                end
                if (vc[3] <= z[1])
                    vc[3] += z[end] - z[1]

                elseif vc[3] >= z[end]
                    vc[3] += z[1] - z[end]
                end
            end
            for i in 1:length(vivecs)
                d = euclid(vc, vivecs[i]);
                if d < min_dist
                    min_dist = d
                    min_idx = i
                end
            end

            if min_dist < 2*dx
                vc = vivecs[min_idx]
                if !vortInBounds(vc, X) && vortInBounds(visorted[end], X)
                    if vc[1] <= x[1]
                        vc[1] += (x[end]-x[1])
                    elseif vc[1] >= x[end]
                        vc[1] += (x[1]-x[end])
                    end
                    if (vc[2] <= y[1])
                        vc[2] += (y[end]-y[1])
    
                    elseif vc[2] >= y[end]
                        vc[2] += (y[1]-y[end])
                    end
                    if (vc[3] <= z[1])
                        vc[3] += z[end] - z[1]
    
                    elseif vc[3] >= z[end]
                        vc[3] += z[1] - z[end]
                    end
                end
                push!(visorted, vc)
                deleteat!(vivecs, min_idx)
            else
                break
            end
        end
        push!(vorts_sorted, visorted)
    end
    return vorts_sorted
end


function sort_classified_vorts2(v_class, vorts_3d, X)
    x = X[1]; y = X[2]; z = X[3];
    dx = x[2]-x[1]; dy = y[2]-y[1]; dz = z[2]-z[1];
    Lx = x[end]-x[1]; Ly = y[end]-y[1]; Lz = z[end]-z[1];
    v_matrix = vcat(vorts_3d'...)[:,1:3]'
    vorts_sorted = []

    for i in 1:length(v_class)
        vi = v_class[i]
        vi = v_matrix[:, collect(vi)]
        vi = vi[:, [vortInBounds2(vi[:, j], X) for j = 1:length(vi[1, :])]] # Filters vortices that aren't on the grid


        vivecs = collect(eachcol(vi))
        visorted_all = []
        visorted = []

        vc = pop!(vivecs)
        v1 = vc
        push!(visorted, vc)
        min_dist = sqrt(3)*(x[end]-x[1]) # init to max distance 
        min_idx = -1

        regimeInBounds = vortInBounds(vc, X)
        changedRegime = false
        while length(vivecs) > 0

            min_dist = 10000
            min_idx = -1
            min_inBounds = false
            ##periodic

            for i in 1:length(vivecs)
                vi_inBounds = vortInBounds(vivecs[i], X)
                if !(regimeInBounds != vi_inBounds && changedRegime)
                    d = euclid(vc, vivecs[i]);
                    if d < min_dist
                        min_dist = d
                        min_idx = i
                        min_inBounds = vi_inBounds
                    end
                end
            end
            jump = false
            if min_dist > 2*dx
                temp = copy(vc)
                if vc[1] <= x[1]
                    temp[1] += Lx
                elseif vc[1] >= x[end]
                    temp[1] -= Lx
                end
                if vc[2] <= y[1]
                    temp[2] += Ly
                elseif vc[2] >= y[end]
                    temp[2] -= Ly
                end
                if vc[3] <= z[1]
                    temp[3] += Lz
                elseif vc[3] >= z[end]
                    temp[3] -= Lz
                end
                for i in 1:length(vivecs)
                    d = euclid(temp, vivecs[i]);
                    if d < min_dist
                        jump = true
                        min_dist = d
                        min_idx = i
                    end
                end
            end 
            if min_dist < 2*dx
                if jump
                    push!(visorted_all, visorted)
                    visorted = []
                    jump  = false
                    regimeInBounds = vortInBounds(vivecs[min_idx], X)
                    changedRegime = false
                
                elseif min_inBounds != regimeInBounds
                    changedRegime = true
                end
                vc = vivecs[min_idx]    
                push!(visorted, vc)
                deleteat!(vivecs, min_idx)
            else
                break
            end
        end
        push!(visorted_all, visorted)
        k = 1
        # while k <= length(visorted_all)
        #     temp = visorted_all[k]
        #     visorted_all[k] = temp[[vortInBounds(temp[j], X) for j = 1:length(temp)]]
        #     if length(visorted_all[k]) == 0
        #         deleteat!(visorted_all, k)
        #         k = k-1
        #     end
        #     k+=1
        # end
        # if euclid(visorted_all[end][end], visorted_all[1][1]) < 3*dx

        #     push!(visorted_all[end], visorted_all[1][1])
        # end
        push!(vorts_sorted, visorted_all)
    end
    return vorts_sorted
end

function sort_classified_vorts3(v_class, vorts_3d, X)
    x = X[1]; y = X[2]; z = X[3];
    dx = x[2]-x[1]; dy = y[2]-y[1]; dz = z[2]-z[1];
    Lx = x[end]-x[1]; Ly = y[end]-y[1]; Lz = z[end]-z[1];
    v_matrix = vcat(vorts_3d'...)[:,1:3]'
    vorts_sorted = []

    for i in 1:length(v_class)
        vi = v_class[i]
        vi = v_matrix[:, collect(vi)]
        vi = vi[:, [vortInBounds3(vi[:, j], X) for j = 1:length(vi[1, :])]] # Filters vortices that aren't on the grid


        vivecs = collect(eachcol(vi))
        visorted = []

        vc = pop!(vivecs)
        push!(visorted, vc)
        min_dist = sqrt(3)*Lx # init to max distance 
        min_idx = -1

        while length(vivecs) > 0

            min_dist = sqrt(3)*(x[end]-x[1]) # init to max distance 
            min_idx = -1
            ##periodic


            for j in 1:length(vivecs)
                v_new = vivecs[j];
                vdx = abs(v_new[1]-vc[1])
                vdy = abs(v_new[2]-vc[2])
                vdz = abs(v_new[3]-vc[3])

                if vdx >= (Lx)
                    vdx -= (Lx+dx)
                end
                if vdy >= (Ly)
                    # println(vc)
                    # println(v_new)
                    vdy -= (Ly+dy)
                end
                if vdz >= (Lz)
                    vdz -= (Lz+dz)
                end

                # d = sqrt(vdx^2 + vdy^2 + vdz^2)
                d = hypot(vdx, vdy, vdz)

                if d < min_dist
                    min_dist = d
                    min_idx = j
                end
            end
            
            if min_dist < dx
                vc = vivecs[min_idx]    
                push!(visorted, vc)
                deleteat!(vivecs, min_idx)
            else
                print("HERE")
                break
            end
        end
        # push!(visorted_all, visorted)

        visorted = visorted[[vortInBounds(visorted[j], X) for j = 1:length(visorted)]]
        push!(vorts_sorted, visorted)
    end
    k = 1
    return vorts_sorted
end
