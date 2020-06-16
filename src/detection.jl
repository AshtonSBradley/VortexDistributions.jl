"""
    vortices = findvortices(ψ<:Field)

Locate vortices as `2π` phase windings around plaquettes on a cartesian spatial grid.

Requires a 2D wavefunction `ψ(x,y)` on a cartesian grid specified by vectors `x`, `y`.

`vortices` - array of vortex coordinates `x,y` and charges `q`.

Each row is of the form `[x, y, cv]`, and the array is sorted into lexical order
according to the vortex `x` coordinates.
"""
function findvortices(psi::Field;periodic=false)
    @unpack ψ,x,y = psi
    vort = findvortices_grid(psi)
    if !periodic
        vort = remove_vortices_edge(vort,psi)
    end

    for (j,vortex) in enumerate(vort)
        v = try
        psi_int,xint,yint = zoom_interp(ψ,x,y,vortex.xv,vortex.yv,periodic=periodic) #periodic: peridic indices here
        v1 = findvortices_grid(Torus(psi_int,xint,yint))
        vint = remove_vortices_edge(v1,Torus(psi_int,xint,yint))[1]
        psi_int,xint,yint = zoom_interp(psi_int,xint,yint,vint.xv,vint.yv)
        v1 = findvortices_grid(Torus(psi_int,xint,yint))
        vint = remove_vortices_edge(v1,Torus(psi_int,xint,yint))[1]
        catch nothing
        end
        v != nothing && (vort[j] = v)    # NOTE fallback to grid if zoom fails
    end
    if periodic
        Lx,Ly = last(x)-first(x),last(y)-first(y)
        vdat = vortex_array(vort)
        xx,yy,qq = vdat[:,1],vdat[:,2],vdat[:,3]
        @. xx[xx>Lx/2] -= Lx
        @. xx[xx<-Lx/2] += Lx
        @. yy[yy>Ly/2] -= Ly
        @. yy[yy<-Lx/2] += Ly
        vort = PointVortex.(xx,yy,qq)
    end
    return vort
end

#TODO add tests for periodic

function findvortices_grid(psi::Torus;shift=true)
    vort = findvortices_jumps(psi,shift=shift)
    varray = vortex_array(vort)
    vort = sortslices(varray,dims=1)
    return PointVortex(vort)
end

function findvortices_grid(psi::Sphere;shift=true)
    @unpack ψ = psi
    windvals = phase_jumps(angle.(ψ),2)
    vort = findvortices_jumps(psi,shift=shift)
    rawvort = vortex_array(vort)

    # North pole winding. - sign means polar vortex co-rotating with w > 0
    w1 = sum(windvals[1,:])
    (sign(w1) != 0) && (rawvort = [rawvort; 0.0 0.0 -w1])
    # South pole winding. + sign means co-rotating with w > 0
    w2 = sum(windvals[end,:])
    (sign(w2) != 0) && (rawvort = [rawvort; pi 0.0 w2])

    vort = sortslices(vort,dims=1)
    return PointVortex(vort)
end

"""
    psiz,xz,yz = zoom(psi,x,y,xv,yv,win=1)

Zoom in near (xv,yv) with a window `win=1`, which gives a 4x4 local grid.
"""
function zoom_grid(psi,x,y,xv,yv;win=1,periodic=false)
    dx,dy = Δ(x),Δ(y)
    nx,ny = size(psi)
    ix = isapprox.(x,xv,atol=dx) |> findfirst
    iy = isapprox.(y,yv,atol=dy) |> findfirst
    ixw = (ix-win):(ix+win+1)
    iyw = (iy-win):(iy+win+1)
    if !periodic
        psiz,xz,yz = psi[ixw,iyw],x[ixw],y[iyw]
    else
        ixp,iyp = mod1.(ixw,nx),mod1.(iyw,ny) # periodic indices
        psiz = psi[ixp,iyp]
        xz = (x[ix]-win*dx):dx:(x[ix]+(win+1)*dx) # local x vector
        yz = (y[iy]-win*dy):dy:(y[iy]+(win+1)*dy) # local y vector
    end
    return psiz,xz,yz
end

"""
    psiz,xz,yz = zoom_interp(psi,x,y,xv,yv,win=1,nz=30)

Zoom in near (xv,yv) with a window `win=1`, which gives a 4x4 local grid.
The wavefunction psi is interpolated onto a 30x30 domain inside the local box.
"""
function zoom_interp(psi,x,y,xv,yv;win=1,nz=30,periodic=false)
    psiz, xz, yz = zoom_grid(psi,x,y,xv,yv,win=win,periodic=periodic)
    psi_itp = interpolate((xz,yz), psiz, Gridded(Linear()))
    xint,yint = LinRange(xz[1],xz[end],nz),LinRange(yz[1],yz[end],nz)
    psi_int = psi_itp(xint,yint)
    return psi_int,xint,yint
end

function findvortices_jumps(psi::Field;shift=true)
    @unpack x,y,ψ = psi
    phase = angle.(ψ)

    # Nearest neighbour phase jumps
    diffx,diffy = phase_jumps(phase,1),phase_jumps(phase,2)

    # Plaquette loop integrals with in-place memory recycling
    circshift!(phase,diffx,(0,1))
    diffx .-= phase; diffx .-= diffy
    circshift!(phase,diffy,(1,0))
    diffx .+= phase

    # Windings
    ixp,iyp,vp = find_where(diffx .> 0.0)
    xp = x[ixp]; yp = y[iyp]; np = length(vp)
    ixn,iyn,vn = find_where(diffx .< 0.0)
    xn = x[ixn]; yn = y[iyn]; nn = length(vn)

    if shift
        dx,dy = Δ(x),Δ(y)
        xp .-= dx/2; yp .-= dy/2; xn .-= dx/2; yn .-= dy/2
    end

    vortices = [xn yn -vn; xp yp vp]

    return PointVortex(vortices)
end
