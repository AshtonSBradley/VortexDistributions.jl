function plot_iso(psi)
    density = abs2.(psi)
    pmax = maximum(density)
    density = density/pmax
    volume(density, algorithm = :iso, show_axis=true)
end

function findvortices3D_z(psi, X)
    x = X[1]; y = X[2];
    z = X[3];
    vorts_3d = []
    for zidx in 1:length(z)
        vorts = vortex_array(findvortices(Torus(psi[:, :, zidx], x, y)));
        if length(vorts) != 0
            for vidx in 1:length(vorts[:, 1])
                v = vorts[vidx, :]
                vz = [v[1], v[2], z[zidx], v[3]]
                push!(vorts_3d, vz)
            end
        end
    end
    return vorts_3d;
end

function findvortices3D_x(psi, X)
    x = X[1]; y = X[2];
    z = X[3];
    vorts_3d = []
    for xidx in 1:length(x)
        vorts = vortex_array(findvortices(Torus(psi[xidx, :, :], y, z)));
        if length(vorts) != 0
            for vidx in 1:length(vorts[:, 1])
                v = vorts[vidx, :]
                vx = [x[xidx], v[1], v[2], v[3]]
                push!(vorts_3d, vx)
            end
        end
    end
    return vorts_3d;
end

function findvortices3D_y(psi, X)
    x = X[1]; y = X[2];
    z = X[3];
    vorts_3d = []
    for yidx in 1:length(y)
        vorts = vortex_array(findvortices(Torus(psi[:, yidx, :], x, z)));
        if length(vorts) != 0
            for vidx in 1:length(vorts[:, 1])
                v = vorts[vidx, :]
                vy = [v[1], y[yidx], v[2], v[3]]
                push!(vorts_3d, vy)
            end
        end
    end
    return vorts_3d;
end

function findidx_uniform(x, arr)
    # For arrays where xi is between [-L, L] and uniform spacing
    # So xi = -L/2 + (L/N)*i
    # => i = (N/L)*(x + L/2)
    N = length(arr);
    L = arr[end]*2;
    idx = (N/L)*(x + L/2);
    if idx < 1
        println("x is not within the bounds of this array");
        return;
    else
        return idx;
    end
end

function findidx_box(x, arr)
    idx = findidx_uniform(x, arr)
    L = 14

    idx = idx/(length(arr)/L);
    idx = idx + 8
    return idx
end 

function plot_vfound3D(vorts_3d, X, m_size=1500)
    x = X[1]; y = X[2]; z = X[3];

    charge = false
    if length(vorts_3d) == 1
        charge = true
    end

    for vorts in vorts_3d
        for vidx in 1:length(vorts)
            v = vorts[vidx]
            vx = findidx_uniform(v[1], x);
            vy = findidx_uniform(v[2], y);
            vz = findidx_uniform(v[3], z);
            vq = v[4]
            if charge
                if vq > 0
                    scatter!([vx], [vy], [vz], color="red", markersize=m_size)
                else
                    scatter!([vx], [vy], [vz], color="blue", markersize=m_size)
                end
            else
                scatter!([vx], [vy], [vz], color="black", markersize=m_size)
            end
        end
    end
end

function plot_vfound3D_box(vorts, X, charge=false)
    x = X[1]; y = X[2]; z = X[3];
    for vidx in 1:length(vorts)
        v = vorts[vidx]
        vx = findidx_uniform(v[1], x);
        vy = findidx_uniform(v[2], y);
        vz = findidx_box(v[3], z);
        vq = v[4]
        if charge
            if vq > 0
                scatter!([vx], [vy], [vz], color="red", markersize=1500)
            else
                scatter!([vx], [vy], [vz], color="blue", markersize=1500)
            end
        else
            scatter!([vx], [vy], [vz], color="black", markersize=1500)
        end
    end
end

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
    return [vorts_xslice, vorts_yslice, vorts_zslice];
end