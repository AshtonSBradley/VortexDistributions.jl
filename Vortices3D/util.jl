
function findvortices3D_itp(psi, X, N=1)
    
    x = X[1]; y = X[2]; z = X[3];

    x_itp = interpolate(X[1], BSpline(Linear()));
    y_itp = interpolate(X[2], BSpline(Linear()));
    z_itp = interpolate(X[3], BSpline(Linear()));

    psi_itp = interpolate(psi, BSpline(Cubic(Line(OnGrid()))))

    x_range = LinRange(1,length(x),N*length(x))
    y_range = LinRange(1,length(y),N*length(y))
    z_range = LinRange(1,length(z),N*length(z))

    vorts_xslice = []
    for xidx in x_range
        vorts = vortex_array(findvortices(Torus(psi_itp[xidx, :, :], y, z)))
        if length(vorts) != 0
            for vidx in 1:length(vorts[:, 1])
                v = vorts[vidx, :]
                vx = [x_itp(xidx), v[1], v[2], v[3]]
                push!(vorts_xslice, vx)
            end
        end
    end

    vorts_yslice = []
    for yidx in y_range
        vorts = vortex_array(findvortices(Torus(psi_itp[:, yidx, :], x, z)))
        if length(vorts) != 0
            for vidx in 1:length(vorts[:, 1])
                v = vorts[vidx, :]
                vy = [v[1], y_itp(yidx), v[2], v[3]]
                push!(vorts_yslice, vy)
            end
        end
    end

    vorts_zslice = []
    for zidx in z_range
        vorts = vortex_array(findvortices(Torus(psi_itp[:, :, zidx], x, y)))
        if length(vorts) != 0
            for vidx in 1:length(vorts[:, 1])
                v = vorts[vidx, :]
                vz = [v[1], v[2], z_itp[zidx], v[3]]
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
    # vorts = vorts3DMatrix(vorts);
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

