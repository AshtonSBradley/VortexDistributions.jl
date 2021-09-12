function find_vortices3D(psi, X)
    x = X[1]; y = X[2];
    Nx = length(x); Ny = length(y);
    if length(X) == 3
        z = X[3];
        Nz = length(z);
 
        for i = 1:Nz
            psi_z = Torus(psi[:,:,i], x, y);

            if i == 1
                vort_array = vortex_array(findvortices(psi_z));
                v = z[i]*ones(length(vort_array[:, 1]));
                vort_array = hcat(vort_array, v);
                vort_array = view(vort_array, :, [1,2,4,3]);
            else 
                temp = vortex_array(findvortices(psi_z));
                v = z[i]*ones(length(temp[:, 1]));
                temp = hcat(temp, v);
                temp = view(temp, :, [1, 2, 4, 3]);
                vort_array = vcat(vort_array, temp);
            end
        end
    else 
        
        vort_array = vortex_array(findvortices(psi));
    end
    temp = zeros(length(vort_array[:, 1]))
    vort_array = hcat(vort_array, temp);
    return vort_array;
end

function find_vortices3D_v2(psi, X)
    x = X[1]; y = X[2];
    z = X[3];
    vorts_3d = Array{Float64, 2}[]
    for zidx in 1:length(z)
        vorts = vortex_array(findvortices(Torus(psi[:, :, zidx], x, y)));
        if (length(vorts) == 0)
            push!(vorts_3d, vorts)
        else
            temp = zeros(length(vorts[:, 1]))
            zval = z[zidx].*ones(length(vorts[:, 1]))
            vorts = hcat(vorts, temp);
            vorts = hcat(vorts, zval);
            vorts =view(vorts, :, [1, 2, 5, 3, 4])
            push!(vorts_3d, vorts);
        end
    end
    return vorts_3d;
end

function find_vortices3D_v2(psi)
    @unpack ψ,x,y,z = psi 
    vorts_3d = Array{Float64, 2}[]
    for zidx in 1:length(z)
        vorts = vortex_array(findvortices(Torus(ψ[:, :, zidx], x, y)));
        if (length(vorts) == 0)
            # push!(vorts_3d, vorts)
        else
            temp = zeros(length(vorts[:, 1]))
            zval = z[zidx].*ones(length(vorts[:, 1]))
            vorts = hcat(vorts, temp);
            vorts = hcat(vorts, zval);
            vorts =view(vorts, :, [1, 2, 5, 3, 4])
            push!(vorts_3d, vorts);
        end
    end
    return vorts_3d;
end
        
# Need to pass X and define x = X[1] etc
function vortex_boundary_bottom(v,X, dx, dy, dz)
    x = X[1]; y = X[2]; z = X[3];
    at_boundary = true;

    if ((v[1, 1] - x[1] > dx) && (x[end] - v[1, 1] > dx)) # Then not at x Boundary
        if ((v[1, 2] - y[1] > dy) && (y[end] - v[1, 2] > dy)) # " y boundary 
            if ((v[1, 3] - z[1] > dz) && (z[end] - v[1, 3] > dz))
                at_boundary = false;
            end
        end
    end
    return at_boundary;
end

function vortex_boundary_top(v, X, dx, dy, dz)
    x = X[1]; y = X[2]; z = X[3];
    at_boundary = true;
    if ((v[end, 1] - x[1] > dx) && (x[end] - v[end, 1] > dx)) # Then not at x Boundary
        if ((v[end, 2] - y[1] > dy) && (y[end] - v[end, 2] > dy)) # " y boundary 
            if ((v[end, 3] - z[1] > dz) && (z[end] - v[end, 3] > dz))
                at_boundary = false;
            end
        end
    end
    return at_boundary;
end

function vortex_link_bottom(vidx, vorts_label, X)
    x = X[1]; y = X[2]; z = X[3];
    dx = x[2]-x[1]; dy = y[2]-y[1]; dz = z[2]-z[1];
    v = vorts_label[vidx];
    for i in 1:length(vorts_label)
        vi = vorts_label[i];
        if vi[1,5] != v[1,5] # different label
            # Find matching z level
            index = 0;
            for j in 1:length(vi[:, 1])
                if vi[j, 3] == v[1, 3]
                    index = j;
                end
            end
            # If matching z level exists
            if index != 0
                # If z level is at bottom of vi 
                if vi[index, 3] == vi[1, 3]
                    vorts_label[i][:, 5] .= vorts_label[vidx][1, 5]
                    if !vortex_boundary_top(vi, X, dx, dy, dz)
                        vortex_link_top(i, vorts_label, X);
                    end
                elseif vi[index, 3] == vi[end, 3]
                    vorts_label[i][:, 5] .= vorts_label[vidx][1, 5]
                    if !vortex_boundary_bottom(vi,X, dx, dy, dz)
                        vortex_link_bottom(i, vorts_label, X);
                    end
                end
            end
        end
    end
end

function vortex_link_top(vidx, vorts_label, X)
    x = X[1]; y = X[2]; z = X[3];
    dx = x[2]-x[1]; dy = y[2]-y[1]; dz = z[2]-z[1];
    v = vorts_label[vidx];
    for i in 1:length(vorts_label)
        vi = vorts_label[i];
        if vi[1,5] != v[1,5] # different label
            # Find matching z level
            index = 0;
            for j in 1:length(vi[:, 1])
                if vi[j, 3] == v[end, 3]
                    index = j;
                end
            end
            # If matching z level exists
            if index != 0
                # If z level is at bottom of vi 
                if vi[index, 3] == vi[1, 3]
                    vorts_label[i][:, 5] .= vorts_label[vidx][1, 5]
                    if !vortex_boundary_top(vi,X, dx, dy, dz)
                        vortex_link_top(i, vorts_label, X);
                    end
                elseif vi[index, 3] == vi[end, 3]
                    vorts_label[i][:, 5] .= vorts_label[vidx][1, 5]

                    if !vortex_boundary_bottom(vi,X, dx, dy, dz)
                        vortex_link_bottom(i, vorts_label, X);
                    end
                end
            end
        end
    end
end

function sort_vorts_label(vorts_label)
    NUM_VORTS = length(vorts_label)
    temp = Array{Float64, 2}[]
    i = 1;
    for i in 1:length(vorts_label)
        for j in 1:NUM_VORTS
        
            if (vorts_label[i][1, 5] == j)
               #println(j)
                push!(temp, vorts_label[i]);
            end
        end
    end
    return temp;
end

function psi_itp_2D(psi_itp, zval, X)
    x = X[1]; y = X[2]; 
    Nx = length(x); Ny = length(y); 
    psi_itp2D = im.*zeros(Nx, Ny);
    for i in 1:Nx
        for j in 1:Ny
                psi_itp2D[i, j] = psi_itp(x[i], y[j], zval);
        end
    end
    return psi_itp2D;   
end