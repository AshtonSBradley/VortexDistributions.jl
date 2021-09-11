using VortexDistributions
using JLD2
using GLMakie
@load "Vortices3D/data/VortexDetection3D/sim_examples/box_rings.jld2" psi_ring1 psi_ring2 psi_ring3 X

## Plot iso surface for select ring 
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
        vorts = vortex_array(findvortices(Torus(psi[xidx, :, :], z, y)));
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

function plot_vfound3D(vorts, X, charge=false)
    x = X[1]; y = X[2]; z = X[3];
    for vidx in 1:length(vorts)
        v = vorts[vidx]
        vx = findidx_uniform(v[1], x);
        vy = findidx_uniform(v[2], y);
        vz = findidx_uniform(v[3], z);
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

psi = psi_ring2;
plot_iso(psi)

vortsx = findvortices3D_x(psi, X);
vortsy = findvortices3D_y(psi, X);
vortsz = findvortices3D_z(psi, X);


plot_vfound3D(vortsx, X)
plot_vfound3D(vortsy, X)
plot_vfound3D(vortsz, X)



