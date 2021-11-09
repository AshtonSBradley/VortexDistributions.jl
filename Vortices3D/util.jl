
function findvortices3D_itp(psi, X, N=1)
    
    x = X[1]; y = X[2]; z = X[3];

    x_itp = interpolate(X[1], BSpline(Linear()));
    y_itp = interpolate(X[2], BSpline(Linear()));
    z_itp = interpolate(X[3], BSpline(Linear()));

    x_etp = extrapolate(x_itp, Line())
    y_etp = extrapolate(y_itp, Line())
    z_etp = extrapolate(z_itp, Line())


    psi_itp = interpolate(psi, BSpline(Quadratic(Periodic(OnCell()))))
    psi_etp = extrapolate(psi_itp, Periodic())

    x_range = LinRange(0,length(x)+1,N*(length(x)))
    y_range = LinRange(0,length(y)+1,N*(length(y)))
    z_range = LinRange(0,length(z)+1,N*(length(z)))

    # x_range = vcat(LinRange(0, 1, 20*N), x_range, LinRange(length(x), length(x)+1, 10*N))
    # y_range = vcat(LinRange(0, 1, 20*N), y_range, LinRange(length(y), length(y)+1, 10*N))
    # z_range = vcat(LinRange(0, 1, 20*N), z_range, LinRange(length(z), length(z)+1, 10*N))

    x = LinRange(-8.25, 8.5, 68);
    y = LinRange(-8.25, 8.5, 68);
    z = LinRange(-8.25, 8.5, 68);

    vorts_xslice = []
    for xidx in x_range
        vorts = vortex_array(findvortices(Torus(psi_etp[xidx, -1:66, -1:66], y, z)))
        if length(vorts) != 0
            for vidx in 1:length(vorts[:, 1])
                v = vorts[vidx, :]
                vx = [x_etp(xidx), v[1], v[2], v[3]]
                push!(vorts_xslice, vx)
            end
        end
    end

    vorts_yslice = []
    for yidx in y_range
        vorts = vortex_array(findvortices(Torus(psi_etp[-1:66, yidx, -1:66], x, z)))
        if length(vorts) != 0
            for vidx in 1:length(vorts[:, 1])
                v = vorts[vidx, :]
                vy = [v[1], y_etp(yidx), v[2], v[3]]
                push!(vorts_yslice, vy)
            end
        end
    end

    vorts_zslice = []
    for zidx in z_range
        vorts = vortex_array(findvortices(Torus(psi_etp[-1:66, -1:66, zidx], x, y)))
        if length(vorts) != 0
            for vidx in 1:length(vorts[:, 1])
                v = vorts[vidx, :]
                vz = [v[1], v[2], z_etp[zidx], v[3]]
                push!(vorts_zslice, vz)
            end
        end
    end
    return vcat([vorts_xslice, vorts_yslice, vorts_zslice]...);
end

function plot_iso(psi, X)
    density = abs2.(psi)
    pmax = maximum(density)
    density = density/pmax
    volume(X[1], X[2], X[3], density, algorithm = :iso, show_axis=true)
end

function scatterVortsOnIso(vorts, markersize=200)
    vorts = vorts3DMatrix(vorts);
    scatter!(vorts[:, 1], vorts[:, 2], vorts[:, 3], color="red", markersize=markersize)
end

function vorts3DMatrix(vorts)
    vorts = vcat(vorts'...)
end

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

function setMethod3(v_matrix, ϵ)
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

function setMethodPeriodic(v_matrix, X, ϵ, periodic=false)
    kdtree = KDTree(v_matrix)

    num_vorts = length(v_matrix[1,:])
    unvisited = Set(collect(1:num_vorts))
    fils = []
    x = X[1]; y = X[2]; z = X[3];
    Δx = x[2]-x[1]; Δy = y[2]-y[1]; Δz = z[2]-z[1];

    xdist = x[end]-x[1]; ydist = y[end]-y[1]; zdist = z[end]-z[1];

    while length(unvisited) > 0
        idx = first(unvisited)
        vc = v_matrix[:, idx]
        f_idxs = inrange(kdtree, vc, ϵ)
        f = Set(f_idxs)
        search = Set(f_idxs)
        setdiff!(search, idx)
        if periodic
            vcx = v_matrix[1,idx]; vcy=v_matrix[2,idx]; vcz = v_matrix[3,idx];
            if abs(vcx - x[1]) < ϵ
                vp = [vcx + xdist + Δx, vcy, vcz]
                p_idxs = inrange(kdtree, vp, ϵ) # Won't include itself this time 
                union!(f, Set(p_idxs))
                union!(search, Set(p_idxs))
            elseif abs(vcx - x[end]) < ϵ
                vp = [vcx - xdist-Δx, vcy, vcz]
                p_idxs = inrange(kdtree, vp, ϵ) # Won't include itself this time 
                union!(f, Set(p_idxs))
                union!(search, Set(p_idxs))
            end
            if abs(vcy - y[1]) < ϵ
                vp = [vcx, vcy + ydist+Δy, vcz]
                p_idxs = inrange(kdtree, vp, ϵ) # Won't include itself this time 
                union!(f, Set(p_idxs))
                union!(search, Set(p_idxs))
            elseif abs(vcy - y[end]) < ϵ
                vp = [vcx, vcy - ydist-Δy, vcz]
                p_idxs = inrange(kdtree, vp, ϵ) # Won't include itself this time 
                union!(f, Set(p_idxs))
                union!(search, Set(p_idxs))
            end
            if abs(vcz - z[1]) < ϵ
                vp = [vcx, vcy, vcz + zdist+Δz]
                p_idxs = inrange(kdtree, vp, ϵ) # Won't include itself this time 
                union!(f, Set(p_idxs))
                union!(search, Set(p_idxs))
            elseif abs(vcz - z[end]) < ϵ
                vp = [vcx, vcy, vcz - zdist-Δz]
                p_idxs = inrange(kdtree, vp,  ϵ) # Won't include itself this time 
                union!(f, Set(p_idxs))
                union!(search, Set(p_idxs))
            end
        end
        while length(search) > 0
            idx = first(search)
            setdiff!(search, idx)
            vc = v_matrix[:, idx]
            vc_idxs = inrange(kdtree, vc, ϵ)
            setdiff!(vc_idxs, f)
            union!(f, Set(vc_idxs))
            union!(search, Set(vc_idxs))
            if periodic
                vcx = v_matrix[1,idx]; vcy=v_matrix[2,idx]; vcz = v_matrix[3,idx];
                if abs(vcx - x[1]) < ϵ
                    vp = [vcx + xdist+Δx, vcy, vcz]
                    p_idxs = inrange(kdtree, vp,  ϵ) # Won't include itself this time
                    setdiff!(p_idxs, f) 
                    union!(f, Set(p_idxs))
                    union!(search, Set(p_idxs))
                elseif abs(vcx - x[end]) < ϵ
                    vp = [vcx - xdist-Δx, vcy, vcz]
                    p_idxs = inrange(kdtree, vp, ϵ) # Won't include itself this time 
                    setdiff!(p_idxs, f) 
                    union!(f, Set(p_idxs))
                    union!(search, Set(p_idxs))
                end
                if abs(vcy - y[1]) < ϵ
                    vp = [vcx, vcy + ydist+Δy, vcz]
                    p_idxs = inrange(kdtree, vp,  ϵ) # Won't include itself this time 
                    setdiff!(p_idxs, f) 
                    union!(f, Set(p_idxs))
                    union!(search, Set(p_idxs))
                elseif abs(vcy - y[end]) < ϵ
                    vp = [vcx, vcy - ydist-Δy, vcz]
                    p_idxs = inrange(kdtree, vp,  ϵ) # Won't include itself this time
                    setdiff!(p_idxs, f) 
                    union!(f, Set(p_idxs))
                    union!(search, Set(p_idxs))
                end
                if abs(vcz - z[1]) < ϵ
                    vp = [vcx, vcy, vcz + zdist+Δz]
                    p_idxs = inrange(kdtree, vp,  ϵ) # Won't include itself this time 
                    setdiff!(p_idxs, f) 
                    union!(f, Set(p_idxs))
                    union!(search, Set(p_idxs))
                elseif abs(vcz - z[end]) < ϵ
                    vp = [vcx, vcy, vcz - zdist-Δz]
                    p_idxs = inrange(kdtree, vp, ϵ) # Won't include itself this time 
                    setdiff!(p_idxs, f) 
                    union!(f, Set(p_idxs))
                    union!(search, Set(p_idxs))
                end
            end
        end
        if length(f) > N
            push!(fils, f)
        end
        setdiff!(unvisited, f)
    end
    return fils
end
            
function naivePlaquette(psi, X, ϵ)
    x = X[1]; y = X[2];
    vorts = []
    phase = angle.(psi)
    for i in 1:(length(x)-1)
        for j in 1:(length(y)-1)
            # start from top left corner and work anti-clockwise
            n = 0
            square = [phase[i, j], phase[i+1, j], phase[i+1, j+1], phase[i, j+1], phase[i, j]]
            for k in 2:5
                diff = square[k]-square[k-1]
                if diff > π
                    count += 1
                    square[k] += count*2*π
                else if diff < - pi
                    count -= 1
                    square[k] += count*2*pi
                end
            end
            phase_diff = square[end] - square[1]
            m = phase_diff/(2*π)
            if abs(m) > ϵ
                ## vort in box
                push!(vorts, [i, j])
            end
        end
    end
end



