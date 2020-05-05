

function findvortices_jumps(psi::Field;shift=true)
    @unpack x,y,ψ = psi
    phase = angle.(ψ)
    dx,dy = x[2]-x[1],y[2]-y[1]
    Lx,Ly = x[end]-x[1]+dx,y[end]-y[1]+dy

    # Compute all nearest neighbour jumps
    diffx = phasejumps(phase,1); diffy = phasejumps(phase,2)

    # Plaquette loop integrals with memory recycling
    circshift!(phase,diffx,(0,1))
    diffx .-= phase; diffx .-= diffy
    circshift!(phase,diffy,(1,0))
    diffx .+= phase

    # Find all windings
    ixp,iyp,vp = findwhere(diffx .> 0.0)
    xp = x[ixp]; yp = y[iyp]; np = length(vp)
    ixn,iyn,vn = findwhere(diffx .< 0.0)
    xn = x[ixn]; yn = y[iyn]; nn = length(vn)

    if shift
        xp .-= dx/2; yp .-= dy/2; xn .-= dx/2; yn .-= dy/2
        @. xp[xp < x[1]] += Lx
        @. xp[xp > x[end]] -= Lx
        @. yp[yp < y[1]] += Ly
        @. yp[yp > y[end]] -= Ly
        @. xn[xn < x[1]] += Lx
        @. xn[xn > x[end]] -= Lx
        @. yn[yn < y[1]] += Ly
        @. yn[yn > y[end]] -= Ly
    end

    vortices = [xn yn -vn; xp yp vp]

    return PointVortex(vortices)
end

function findvortices_grid(psi::Torus;shift=true)
    vort = findvortices_jumps(psi,shift=shift)
    rawvort = rawData(vort)
    vort = sortslices(rawvort,dims=1)
    return PointVortex(vort)
end

function findvortices_grid(psi::Sphere;shift=true)
    vort = findvortices_jumps(psi,shift=shift)
    rawvort = rawData(vort)

    @unpack ψ = psi
    # North pole winding. - sign means polar vortex co-rotating with w > 0
    windvals = phasejumps(angle.(ψ[1,:]))
    w1 = sum(windvals)
    (sign(w1) != 0) && (rawvort = [rawvort; 0.0 0.0 -w1])
    # South pole winding. + sign means co-rotating with w > 0
    windvals = phasejumps(angle.(ψ[end,:]))
    w2 = sum(windvals)
    (sign(w2) != 0) && (rawvort = [rawvort; pi 0.0 w2])

    vort = sortslices(vort,dims=1)
    return PointVortex(vort)
end

function findvortices_interp(psi::Field)
    vort = findvortices_grid(psi,shift=true)
    # vort = removeedgevortices(vort,psi) # TODO periodic detection
    dx = psi.x[2]-psi.x[1]
    dy = psi.y[2]-psi.y[1]
    Lx = psi.x[end]-psi.x[1]+dx
    Ly = psi.y[end]-psi.y[1]+dy

# TODO: replace corezoom with another sub-grid method
    for (j,vortex) in enumerate(vort)
        try
        vortz,psiz = corezoom(vortex,psi)
        vortz,psiz = corezoom(vortz,psiz)
        vortz,psiz = corezoom(vortz,psiz)
        vort[j] = vortz
        catch nothing
        end
        # put coordinates in periodic form
        vort[j].xv < psi.x[1] && (vort[j].xv += Lx)
        vort[j].xv > psi.x[end] && (vort[j].xv -= Lx)
        vort[j].yv < psi.y[1] && (vort[j].yv += Ly)
        vort[j].yv > psi.y[end] && (vort[j].yv -= Ly)
    end

    return vort
end

"""
    vortices = findvortices(ψ<:Field;interp=true)

Locate vortices as `2π` phase windings around plaquettes on a cartesian spatial grid.

Requires a 2D wavefunction `ψ(x,y)` on a cartesian grid specified by vectors `x`, `y`.

`vortices` - array of vortex coordinates `xv,yv` and charges `qv`.

Each row is of the form `[xv, yv, cv]`, and the array is sorted into lexical order
according to the `xv` coordinates.
"""
function findvortices(psi::Field;interp=true)
    if interp
        return findvortices_interp(psi)
    else
        return findvortices_grid(psi)
    end
end

"""
    vortz,psiz = corezoom(vortex,ψ<:Field,winhalf=2,Nz=30)

Uses local interpolation to resolve core location.
"""
function corezoom(vortex::PointVortex,psi::T,winhalf=2,Nz=30) where T<:Field
    @unpack ψ,x,y = psi
    s1,s2 = size(psi)
    xv,yv,qv = rawData(vortex)
    dx=x[2]-x[1]; dy=y[2]-y[1]
    ixv = isapprox.(x,xv,atol=dx) |> findlast
    iyv = isapprox.(y,yv,atol=dy) |> findlast
    ixwin = (ixv-winhalf):(ixv+winhalf-1)
    iywin = (iyv-winhalf):(iyv+winhalf-1)
    xw = (x[ixv]-winhalf*dx):(x[ixv]+(winhalf-1)*dx)
    yw = (y[iyv]-winhalf*dy):(y[iyv]+(winhalf-1)*dy)
    psiw = ψ[mod1(ixwin,end),mod1(iywin,end)]
    xz = LinRange(xw[1],xw[end],Nz)
    yz = LinRange(yw[1],yw[end],Nz)
    knots = (xw,yw)
    itp = interpolate(knots, psiw, Gridded(Linear()))
    psiz = itp(xz,yz)
    ψv = T(psiz,xz |> Vector,yz |> Vector)
    vortz = findvortices_grid(ψv,shift=false)
    # vortz = removeedgevortices(vortz,ψv)[1]
    return vortz,ψv
end

corezoom(vortex::Array{PointVortex,1},psi::Field,winhalf=2,Nz=30) = corezoom(vortex[1],psi,2,30)
