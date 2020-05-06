"""
    vortices = findvortices(ψ<:Field,interp=true)

Locate vortices as `2π` phase windings around plaquettes on a cartesian spatial grid.

Requires a 2D wavefunction `ψ(x,y)` on a cartesian grid specified by vectors `x`, `y`.

`vortices` - array of vortex coordinates `xv,yv` and charges `qv`.

Each row is of the form `[xv, yv, cv]`, and the array is sorted into lexical order
according to the `xv` coordinates.
"""
function findvortices(psi::Field,interp=true)
    if interp
        return findvortices_interp(psi)
    else
        return findvortices_grid(psi)
    end
end

function findvortices_grid(psi::Torus;shift=true)
    vort = findvortices_jumps(psi,shift=shift)
    rawvort = rawData(vort)
    vort = sortslices(rawvort,dims=1)
    return PointVortex(vort)
end

function findvortices_grid(psi::Sphere;shift=true)
    @unpack ψ = psi
    windvals = phasejumps(angle.(ψ),2)
    vort = findvortices_jumps(psi,shift=shift)
    rawvort = rawData(vort)

    # North pole winding. - sign means polar vortex co-rotating with w > 0
    w1 = sum(windvals[1,:])
    (sign(w1) != 0) && (rawvort = [rawvort; 0.0 0.0 -w1])
    # South pole winding. + sign means co-rotating with w > 0
    w2 = sum(windvals[end,:])
    (sign(w2) != 0) && (rawvort = [rawvort; pi 0.0 w2])

    vort = sortslices(vort,dims=1)
    return PointVortex(vort)
end

function findvortices_interp(psi::Field)
    vort = findvortices_grid(psi,shift=true)
    vort = removeedgevortices(vort,psi) # essential for current foundNear!
    #TODO: allow for interp with periodic data (here edges are stripped)
    #NOTE: foundNear uses randVortexField which requires removing endgevortices
    #to make periodic we have to remove all statements that remove edge vortices,
    #except inside ingerpolation which will generate spurious edge vortices.

    for (j,vortex) in enumerate(vort)
        try
        vortz,psiz = corezoom(vortex,psi)
        # vortz,psiz = corezoom_periodic(vortex,psi)
        vortz,psiz = corezoom(vortz,psiz)
        vortz,psiz = corezoom(vortz,psiz)
        vort[j] = vortz
        catch nothing
        end
    end
    return vort
end

function findvortices_jumps(psi::Field;shift=true)
    @unpack x,y,ψ = psi
    phase = angle.(ψ)

    # Nearest neighbour phase jumps
    diffx = phasejumps(phase,1); diffy = phasejumps(phase,2)

    # Plaquette loop integrals with in-place memory recycling
    circshift!(phase,diffx,(0,1))
    diffx .-= phase; diffx .-= diffy
    circshift!(phase,diffy,(1,0))
    diffx .+= phase

    # Windings
    ixp,iyp,vp = findwhere(diffx .> 0.0)
    xp = x[ixp]; yp = y[iyp]; np = length(vp)
    ixn,iyn,vn = findwhere(diffx .< 0.0)
    xn = x[ixn]; yn = y[iyn]; nn = length(vn)

    if shift
        dx,dy=Δ(x),Δ(y)
        xp .-= dx/2; yp .-= dy/2; xn .-= dx/2; yn .-= dy/2
        Lx = x[end]-x[1]+dx
        Ly = y[end]-y[1]+dy
        @. xp[xp>Lx/2] -= Lx
        @. xp[xp<-Lx/2] += Lx
        @. xn[xn>Lx/2] -= Lx
        @. xn[xn<-Lx/2] += Lx
        @. yp[yp>Ly/2] -= Ly
        @. yp[yp<-Ly/2] += Ly
        @. yn[yn>Ly/2] -= Ly
        @. yn[yn<-Ly/2] += Ly
    end

    vortices = [xn yn -vn; xp yp vp]

    return PointVortex(vortices)
end


"""
    vortz,psiz = corezoom(vortex,ψ<:Field,winhalf=2,Nz=30)

Uses local interpolation to resolve core location.
"""
function corezoom(vortex::PointVortex,psi::T,winhalf=2,Nz=30) where T<:Field
    @unpack ψ,x,y = psi
    xv,yv,qv = rawData(vortex)
    dx,dy=Δ(x),Δ(y)
    ixv = isapprox.(x,xv,atol=dx) |> findfirst
    iyv = isapprox.(y,yv,atol=dy) |> findfirst
    ixwin = (ixv-winhalf):(ixv+winhalf-1)
    iywin = (iyv-winhalf):(iyv+winhalf-1)
    xw = x[ixwin]; yw = y[iywin]; psiw = ψ[ixwin,iywin]
    xz = LinRange(xw[1],xw[end],Nz)
    yz = LinRange(yw[1],yw[end],Nz)
    knots = (xw,yw)
    itp = interpolate(knots, psiw, Gridded(Linear()))
    psiz = itp(xz,yz)
    ψv = T(psiz,xz |> Vector,yz |> Vector)
    vortz = findvortices_grid(ψv,shift=false)
    vortz = removeedgevortices(vortz,ψv)[1]
    return vortz,ψv
end

corezoom(vortex::Array{PointVortex,1},psi::Field,winhalf=2,Nz=30) = corezoom(vortex[1],psi,2,30)

function corezoom_periodic(vortex::PointVortex,psi::T,winhalf=2,Nz=30) where T<:Field
    @unpack ψ,x,y = psi
    xv,yv,qv = rawData(vortex)
    dx,dy = Δ(x),Δ(y)
    ixv = isapprox.(x,xv,atol=dx) |> findfirst
    iyv = isapprox.(y,yv,atol=dy) |> findfirst
    ixwin = (ixv-winhalf):(ixv+winhalf-1)
    iywin = (iyv-winhalf):(iyv+winhalf-1)
    xw = (x[ixv]-winhalf*dx):(x[ixv]+winhalf*dx-dx)
    yw = (y[iyv]-winhalf*dy):(y[iyv]+winhalf*dy-dy)
    psiw = ψ[mod1.(ixwin,end),mod1.(iywin,end)]
    knots = (xw,yw)
    itp = interpolate(knots, psiw, Gridded(Linear()))
    xz = LinRange(xw[1],xw[end],Nz)
    yz = LinRange(yw[1],yw[end],Nz)
    psiz = itp(xz,yz)
    ψv = T(psiz,xz |> Vector,yz |> Vector)
    vortz = findvortices_grid(ψv,shift=false)
    vortz = removeedgevortices(vortz,ψv)[1]
    return vortz,ψv
end

corezoom_periodic(vortex::Array{PointVortex,1},psi::Field,winhalf=2,Nz=30) = corezoom_periodic(vortex[1],psi,2,30)
