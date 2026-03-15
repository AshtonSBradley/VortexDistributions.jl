function findvortices(psi::Field; periodic = false, backend::AbstractDetectionBackend2D = default_backend2d(psi))
    @unpack ψ, x, y = psi
    vort = findvortices_grid(psi; backend = backend)
    if !periodic
        vort = remove_vortices_edge(vort, psi)
    end

    for (j, vortex) in enumerate(vort)
        v = try
            psi_int, xint, yint = zoom_interp(ψ, x, y, vortex.xv, vortex.yv, periodic = periodic)
            v1 = findvortices_grid(Torus(psi_int, xint, yint); backend = backend)
            vint = remove_vortices_edge(v1, Torus(psi_int, xint, yint))[1]
            psi_int, xint, yint = zoom_interp(psi_int, xint, yint, vint.xv, vint.yv)
            v1 = findvortices_grid(Torus(psi_int, xint, yint); backend = backend)
            remove_vortices_edge(v1, Torus(psi_int, xint, yint))[1]
        catch
            nothing
        end
        v !== nothing && (vort[j] = v)
    end

    if periodic
        Lx, Ly = last(x) - first(x), last(y) - first(y)
        vdat = vortex_array(vort)
        if size(vdat) != (0, 0)
            xx, yy, qq = vdat[:, 1], vdat[:, 2], vdat[:, 3]
            @. xx[xx > Lx / 2] -= Lx
            @. xx[xx < -Lx / 2] += Lx
            @. yy[yy > Ly / 2] -= Ly
            @. yy[yy < -Ly / 2] += Ly
            vort = PointVortex.(xx, yy, qq)
        end
    end
    return vort
end

function findvortices_grid(psi::Field; shift = true, backend::AbstractDetectionBackend2D = default_backend2d(psi))
    return _findvortices_grid(backend, psi; shift = shift)
end

function _findvortices_grid(::PlaquetteBackend2D, psi::Torus; shift = true)
    vort = findvortices_jumps(psi; shift = shift)
    varray = vortex_array(vort)
    vort = sortslices(varray, dims = 1)
    return PointVortex(vort)
end

function _findvortices_grid(::PlaquetteBackend2D, psi::Sphere; shift = true)
    @unpack ψ = psi
    windvals = phase_jumps(angle.(ψ), 2)
    vort = findvortices_jumps(psi; shift = shift)
    rawvort = vortex_array(vort)

    w1 = sum(windvals[1, :])
    sign(w1) != 0 && (rawvort = [rawvort; 0.0 0.0 -w1])

    w2 = sum(windvals[end, :])
    sign(w2) != 0 && (rawvort = [rawvort; pi 0.0 w2])

    vort = sortslices(rawvort, dims = 1)
    return PointVortex(vort)
end

function findvortices_jumps(psi::Field; shift = true, backend::AbstractDetectionBackend2D = default_backend2d(psi))
    return _findvortices_jumps(backend, psi; shift = shift)
end

function _findvortices_jumps(::PlaquetteBackend2D, psi::Field; shift = true)
    @unpack x, y, ψ = psi
    phase = angle.(ψ)

    diffx, diffy = phase_jumps(phase, 1), phase_jumps(phase, 2)

    circshift!(phase, diffx, (0, 1))
    diffx .-= phase
    diffx .-= diffy
    circshift!(phase, diffy, (1, 0))
    diffx .+= phase

    ixp, iyp, vp = find_where(diffx .> 0.0)
    xp = x[ixp]
    yp = y[iyp]
    ixn, iyn, vn = find_where(diffx .< 0.0)
    xn = x[ixn]
    yn = y[iyn]

    if shift
        dx, dy = Δ(x), Δ(y)
        xp .-= dx / 2
        yp .-= dy / 2
        xn .-= dx / 2
        yn .-= dy / 2
    end

    vortices = [xn yn -vn; xp yp vp]
    return PointVortex(vortices)
end

function found_near(n)
    near = true
    for _ in 1:n
        psi, vort = rand_vortexfield(1)
        vortfound = findvortices(psi)
        vfdata = vortex_array(vortfound)
        vdata = vortex_array(vort)
        dx = Δ(psi.x)
        near *= isapprox(vdata, vfdata, rtol = dx / 4)
    end
    return near
end

function found_near(n, nvorts)
    near = true
    for _ in 1:n
        psi, vort = rand_vortexfield(nvorts)
        vortfound = findvortices(psi)
        vfdata = vortex_array(vortfound)
        vdata = vortex_array(vort)
        if length(vfdata[:, 1]) != nvorts
            return false
        end
        dx = Δ(psi.x)
        sort!(vdata, dims = 1)
        sort!(vfdata, dims = 1)
        for i in 1:nvorts
            near *= isapprox(vdata[i, :], vfdata[i, :], rtol = dx / 4)
        end
    end
    return near
end

function remove_vortices_edge(vort::Array{PointVortex, 1}, psi::Field, edge = 1)
    @unpack x, y = psi
    dx, dy = Δ(x), Δ(y)
    keep = Int[]
    for j in 1:length(vort)
        xi, yi, _ = vortex_array(vort[j])
        xedge = isapprox(xi, x[1], atol = edge * dx) || isapprox(xi, x[end], atol = edge * dx)
        yedge = isapprox(yi, y[1], atol = edge * dy) || isapprox(yi, y[end], atol = edge * dy)
        !(xedge || yedge) && push!(keep, j)
    end
    return vort[keep]
end
