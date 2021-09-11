# Compute the gradient of 2D array, returns grad = [real, imag]

function gradient_2D(arr, X)
    arr_real = real(arr)
    arr_imag = imag(arr)
    x = X[1]; y = X[2];
    dx = x[2] - x[1];
    dy = y[2] - y[1];
    dxarr_real = similar(arr_real);
    dyarr_real = similar(arr_real); 
    dxarr_imag = similar(arr_imag);
    dyarr_imag = similar(arr_imag);

    for xin in 2:length(x)
        dxarr_real[xin, :] .= (arr_real[xin, :] .- arr_real[xin-1, :]) ./ dx;
        dxarr_imag[xin, :] .= (arr_imag[xin, :] .- arr_imag[xin-1, :]) ./ dx;
    end
    dxarr_real[1, :] .= 0;
    dxarr_imag[1, :] .= 0;
    
    for yin in 2:length(y)
        dyarr_real[:, yin] .= (arr_real[:, yin] .- arr_real[:, yin-1]) ./ dy;
        dyarr_imag[:, yin] .= (arr_imag[:, yin] .- arr_imag[:, yin-1]) ./ dy;
    end
    dyarr_real[:, 1] .= 0;
    dyarr_imag[:, 1] .= 0;
    GRADarr_real = [dxarr_real, dyarr_real];
    GRADarr_imag = [dxarr_imag, dyarr_imag];
    grad = [GRADarr_real, GRADarr_imag];
    return grad;
end

# Compute the gradient of 3D array, returns grad = [real, imag]
function gradient_3D(arr, X)

    arr_real = real(arr)
    arr_imag = imag(arr)
    x = X[1]; y = X[2]; z = X[3];
    dx = x[2] - x[1];
    dy = y[2] - y[1];
    dz = z[2] - z[1];
  
    dxarr_real = similar(arr_real);
    dyarr_real = similar(arr_real); 
    dzarr_real = similar(arr_real);
    dxarr_imag = similar(arr_imag);
    dyarr_imag = similar(arr_imag);
    dzarr_imag = similar(arr_imag); 

    for xin in 2:length(x)
        dxarr_real[xin, :, :] .= (arr_real[xin, :, :] .- arr_real[xin-1, :, :]) ./ dx;
        dxarr_imag[xin, :, :] .= (arr_imag[xin, :, :] .- arr_imag[xin-1, :, :]) ./ dx;
    end
    dxarr_real[1, :, :] .= 0;
    dxarr_imag[1, :, :] .= 0;

    for yin in 2:length(y)
        dyarr_real[:, yin, :] .= (arr_real[:, yin, :] .- arr_real[:, yin-1, :]) ./ dy;
        dyarr_imag[:, yin, :] .= (arr_imag[:, yin, :] .- arr_imag[:, yin-1, :]) ./ dy;
    end
    dyarr_real[:, 1, :] .= 0;
    dyarr_imag[:, 1, :] .= 0;

    for zin in 2:length(z)
        dzarr_real[:, :, zin] .= (arr_real[:, :, zin] .- arr_real[:, :, zin - 1]) ./ dz;
        dzarr_imag[:, :, zin] .= (arr_imag[:, :, zin] .- arr_imag[:, :, zin - 1]) ./ dz;
    end
    dzarr_real[:, :, 1] .= 0;
    dzarr_imag[:, :, 1] .= 0;

    GRADarr_real = [dxarr_real, dyarr_real, dzarr_real];
    GRADarr_imag = [dxarr_imag, dyarr_imag, dzarr_imag];
    GRADpsi = [GRADarr_real, GRADarr_imag];
    return GRADpsi;
end

# Pseudo-vorticity 3D, needs X = [x, y, z] for the grid and array gradient as parameters.
# Returns pseudo-vorticity vector 
function wps_v(gradpsi, X)
    x = X[1]; y = X[2]; z = X[3];
    grad_real = gradpsi[1]; grad_imag = gradpsi[2];
    dx_r = grad_real[1]; dy_r = grad_real[2]; dz_r = grad_real[3];
    dx_i = grad_imag[1]; dy_i = grad_imag[2]; dz_i = grad_imag[3];
    grad_v_r = [dx_r[x,y,z], dy_r[x,y,z], dz_r[x,y,z]];
    grad_v_i = [dx_i[x,y,z], dy_i[x,y,z], dz_i[x,y,z]];


    cross_v = cross(grad_v_r,  grad_v_i)
    wps = cross_v
    return wps;
end

# Pseudo-vorticity 2D, needs X = [x, y] for the grid and array gradient as parameters.
# Returns vector orthog to x, y. Good to sanity check as wps = 0 if no nearby vorticies.
function wps_v_2D(gradpsi, X)
    x = X[1]; y = X[2];
    grad_real = gradpsi[1]; grad_imag = gradpsi[2];
    dx_r = grad_real[1]; dy_r = grad_real[2];
    dx_i = grad_imag[1]; dy_i = grad_imag[2];
    grad_v_r = [dx_r[x,y], dy_r[x,y], 0];
    grad_v_i = [dx_i[x,y], dy_i[x,y], 0];

    cross_v = cross(grad_v_r, grad_v_i)
    wps = cross_v;
    return wps;
end


# Find closest element to x in array
function find_nearest(x, arr)
    nearest = arr[1];
    nearest_idx = 1;
    diff = abs(x-nearest);

    for i in 2:length(arr)
        if abs(x-arr[i]) < diff
            nearest = arr[i];
            nearest_idx = i;
            diff = abs(x-arr[i]);
        end
    end
    return [nearest,  Int(nearest_idx)];
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

function wps_array_plot(psi, zslice)
    grad = gradient_3D_cent(psi, X);
    psi_torus = Torus(psi[:, :, zslice], x, y);
    vorts = vortex_array(findvortices(psi_torus))
    wps_arr = similar(real(psi[:, :, 1]))

    for i in 1:length(x)
        for j in 1:length(y)
            wps_arr[i, j] = norm(wps_v(grad, [i,j, zslice]))
        end
    end
    wps_arr = abs.(wps_arr)
    max = findmax(wps_arr)[1]
    wps_arr = wps_arr ./max

    heatmap(x, y, wps_arr)
    scatter!(vorts[:, 2], vorts[:, 1])
end

function gradient_3D_cent(arr, X)
    x = X[1]; y = X[2]; z = X[3];
    arr_r = real(arr);
    arr_i = imag(arr);

    dx = x[2] - x[1]; dy = y[2] - y[1]; dz = z[2] - z[1];
    dxarr_r = similar(arr_r);
    dyarr_r = similar(arr_r); 
    dzarr_r = similar(arr_r);
    dxarr_i = similar(arr_i);
    dyarr_i = similar(arr_i);
    dzarr_i = similar(arr_i); 

    # Periodic Boundary conditions
    dxarr_r[1, :, :] .= (arr_r[2, :, :] .- arr_r[end, :, :]) ./ (2*dx);
    dxarr_r[end, :, :] .= (arr_r[1, :, :] .- arr_r[end-1, :, :]) ./ (2*dx);
    dxarr_i[1, :, :] .= (arr_i[2, :, :] .- arr_i[end, :, :]) ./ (2*dx);
    dxarr_i[end, :, :] .= (arr_i[1, :, :] .- arr_i[end-1, :, :]) ./ (2*dx);

    for i in 2:length(x)-1
        dxarr_r[i,:,:] .= (arr_r[i+1, :, :] .- arr_r[i-1, :, :]) ./ (2*dx);
        dxarr_i[i,:,:] .= (arr_i[i+1, :, :] .- arr_i[i-1, :, :]) ./ (2*dx);
    end

    dyarr_r[:, 1, :] .= (arr_r[:, 2, :] .- arr_r[:, end, :]) ./ (2*dy);
    dyarr_r[:, end, :] .= (arr_r[:, 1, :] .- arr_r[:, end-1, :]) ./ (2*dy);
    dyarr_i[:, 1, :] .= (arr_i[:, 2, :] .- arr_i[:, end, :]) ./ (2*dy);
    dyarr_i[:, end, :] .= (arr_i[:, 1, :] .- arr_i[:, end-1, :]) ./ (2*dy);

    for i in 2:length(y)-1
        dyarr_r[:,i,:] .= (arr_r[:, i+1, :] .- arr_r[:, i-1, :]) ./ (2*dy);
        dyarr_i[:,i,:] .= (arr_i[:, i+1, :] .- arr_i[:, i-1, :]) ./ (2*dy);
    end

    dzarr_r[:, :, 1] .= (arr_r[:, :, 2] .- arr_r[:, :, end])./ (2*dz);
    dzarr_r[:, :, end] .= (arr_r[:, :, 1] .- arr_r[:, :, end-1])./ (2*dz);
    dzarr_i[:, :, 1] .= (arr_i[:, :, 2] .- arr_i[:, :, end])./ (2*dz);
    dzarr_i[:, :, end] .= (arr_i[:, :, 1] .- arr_i[:, :, end-1])./ (2*dz);

    for i in 2:length(z)-1
        dzarr_r[:, :, i] .= (arr_r[:, :, i+1] - arr_r[:, :, i-1]) ./ (dz*2);
        dzarr_i[:, :, i] .= (arr_i[:, :, i+1] - arr_i[:, :, i-1]) ./ (dz*2);
    end
    
    grad_r = [dxarr_r, dyarr_r, dzarr_r];
    grad_i = [dxarr_i, dyarr_i, dzarr_i];
    grad = [grad_r, grad_i];

    return grad;
end

function grad_itp(grad, X)
    x = X[1]; y = X[2]; z = X[3];
    x = LinRange(x[1], x[end], length(x));
    y = LinRange(y[1], y[end], length(y));
    z = LinRange(z[1], z[end], length(z));
    grad_r = grad[1]; grad_i = grad[2];
    grx = grad_r[1];
    gry = grad_r[2];
    grz = grad_r[3];
    gix = grad_i[1];
    giy = grad_i[2];
    giz = grad_i[3];
    grx_itp = interpolate(grx, BSpline(Cubic(Line(OnGrid()))));
    gry_itp = interpolate(gry, BSpline(Cubic(Line(OnGrid()))));
    grz_itp = interpolate(grz, BSpline(Cubic(Line(OnGrid()))));
    gix_itp = interpolate(gix, BSpline(Cubic(Line(OnGrid()))));
    giy_itp = interpolate(giy, BSpline(Cubic(Line(OnGrid()))));
    giz_itp = interpolate(giz, BSpline(Cubic(Line(OnGrid()))));
    sgrx_itp = scale(grx_itp, x, y, z);
    sgry_itp = scale(gry_itp, x, y, z);
    sgrz_itp = scale(grz_itp, x, y, z);
    sgix_itp = scale(gix_itp, x, y, z);
    sgiy_itp = scale(giy_itp, x, y, z);
    sgiz_itp = scale(giz_itp, x, y, z);
    gr_v = [sgrx_itp, sgry_itp, sgrz_itp];
    gi_v = [sgix_itp, sgiy_itp, sgiz_itp];
    grad = [gr_v, gi_v];
    return grad;
end

function wps_int(grad, Xin)
   
    grad_r = grad[1]; grad_i = grad[2];

    gr_v = [grad_r[1](Xin[1], Xin[2], Xin[3]), grad_r[2](Xin[1], Xin[2], Xin[3]), grad_r[3](Xin[1], Xin[2], Xin[3])];
    gi_v = [grad_i[1](Xin[1], Xin[2], Xin[3]), grad_i[2](Xin[1], Xin[2], Xin[3]),  grad_i[3](Xin[1], Xin[2], Xin[3])];

    wps = cross(gr_v, gi_v);
    return wps;
end

function psi_itp(psi, X)
    x = X[1]; y = X[2]; z = X[3];
    x = LinRange(x[1], x[end], length(x));
    y = LinRange(y[1], y[end], length(y));
    z = LinRange(z[1], z[end], length(z));
    psi_itp = interpolate(psi, BSpline(Cubic(Line(OnGrid()))));
    psi_itp = scale(psi_itp, x, y, z);

    return psi_itp
end

function vortices_plot(psi, x, y, zidx)
    psi = Torus(psi_xspc[:, :, zidx], x, y);
    vort_array = vortex_array(findvortices(psi));
    vort_pos = vort_array[:, 1:2];
 
    @views Plots.heatmap(x, y, angle.(psi_xspc[:, :, zidx]), title = "z level = $zidx", label = "Phase")
    Plots.scatter!(vort_pos[:, 2], vort_pos[:, 1], color =:white, label = "Vortices")
    Plots.xlabel!("x")
    Plots.ylabel!("y")

    #when calling:
    #=
    gr()
    anim = @animate for zidx = 1:Nz
        vortices_plot(psi_kspc, x, y, zidx)
    end
    gif(anim, "anim_fps15.gif", fps = 8)
    =#

end

function find_closest_tuple(arr, X)
    x = X[1]; y = X[2];
    closest_dist = 1000;
    index = 0;
    for i in 1:length(arr[:, 1])
        dist = (x - arr[i, 1])^2 + (y - arr[i, 2])^2;
        dist = sqrt(dist);
        if dist < closest_dist
            index = i;
            closest_dist = dist;
        end
    end

    return index;
end