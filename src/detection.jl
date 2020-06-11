"""
    vortices = findvortices(ψ<:Field)

Locate vortices as `2π` phase windings around plaquettes on a cartesian spatial grid.

Requires a 2D wavefunction `ψ(x,y)` on a cartesian grid specified by vectors `x`, `y`.

`vortices` - array of vortex coordinates `x,y` and charges `q`.

Each row is of the form `[x, y, cv]`, and the array is sorted into lexical order
according to the vortex `x` coordinates.
"""
function findvortices(psi::Field)
    @unpack ψ,x,y = psi
    vort = findvortices_grid(psi)
    vort = remove_vortices_edge(vort,psi)

    for (j,vortex) in enumerate(vort)
        v = try
        psi_int,xint,yint = zoom_interp(ψ,x,y,vortex.xv,vortex.yv)
        v1 = findvortices_grid(Torus(psi_int,xint,yint))
        vint = remove_vortices_edge(v1,Torus(psi_int,xint,yint))[1]
        psi_int,xint,yint = zoom_interp(psi_int,xint,yint,vint.xv,vint.yv)
        v1 = findvortices_grid(Torus(psi_int,xint,yint))
        vint = remove_vortices_edge(v1,Torus(psi_int,xint,yint))[1]
        catch nothing
        end
        vort[j] = v     # TODO a fallback for not found
    end
    return vort
end

#TODO essential for current tests: found_near!
#TODO: allow for interp with periodic data (here edges are stripped)
#NOTE: found_near uses rand_vortexfield which requires removing endge vortices
#to make periodic we have to remove all statements that remove edge vortices,
#except inside interpolation which will generate spurious edge vortices.

# first zoom
# psi_int,xint,yint = zoom_interp(psi.ψ,x,y,vortex_array(vp)[1:2]...)
# v1 = findvortices_grid(Torus(psi_int,xint,yint))
# vint = remove_vortices_edge(v1,Torus(psi_int,xint,yint))[1]
#
# # second zoom
# psi_int,xint,yint = zoom_interp(psi_int,xint,yint,vint.xv,vint.yv)
# v2 = findvortices_grid(Torus(psi_int,xint,yint))
# vint = remove_vortices_edge(v2,Torus(psi_int,xint,yint))[1]

# function findvortices(psi::Field)
#     vort = findvortices_grid(psi,shift=true)
#     vort = remove_vortices_edge(vort,psi)
#     #TODO essential for current tests: found_near!
#     #TODO: allow for interp with periodic data (here edges are stripped)
#     #NOTE: found_near uses rand_vortexfield which requires removing endge vortices
#     #to make periodic we have to remove all statements that remove edge vortices,
#     #except inside interpolation which will generate spurious edge vortices.
#
#     for (j,vortex) in enumerate(vort)
#         try
#         vortz,psiz = corezoom(vortex,psi)
#         # vortz,psiz = corezoom_periodic(vortex,psi) #TODO
#         vortz,psiz = corezoom(vortz,psiz)
#         vortz,psiz = corezoom(vortz,psiz)
#         vort[j] = vortz
#         catch nothing
#         end
#     end
#     return vort
# end

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
        #TODO handle periodicity
        # Lx = x[end]-x[1]
        # Ly = y[end]-y[1]
        # @. xp[xp>Lx/2] -= Lx
        # @. xp[xp<-Lx/2] += Lx
        # @. xn[xn>Lx/2] -= Lx
        # @. xn[xn<-Lx/2] += Lx
        # @. yp[yp>Ly/2] -= Ly
        # @. yp[yp<-Ly/2] += Ly
        # @. yn[yn>Ly/2] -= Ly
        # @. yn[yn<-Ly/2] += Ly
    end

    vortices = [xn yn -vn; xp yp vp]

    return PointVortex(vortices)
end

"""
    psiz,xz,yz = zoom(psi,x,y,xv,yv,win=1)

Zoom in near (xv,yv) with a window `win=1`, which gives a 4x4 local grid.
"""
function zoom_grid(psi,x,y,xv,yv;win=1)
    dx,dy = Δ(x),Δ(y)
    ix = isapprox.(x,xv,atol=dx) |> findfirst
    iy = isapprox.(y,yv,atol=dy) |> findfirst
    ixw = (ix-win):(ix+win+1)
    iyw = (iy-win):(iy+win+1)
    psiz,xz,yz = psi[ixw,iyw],x[ixw],y[iyw]
    return psiz,xz,yz
end

"""
    psiz,xz,yz = zoom_interp(psi,x,y,xv,yv,win=1,nz=30)

Zoom in near (xv,yv) with a window `win=1`, which gives a 4x4 local grid.
The wavefunction psi is interpolated onto a 30x30 domain inside the local box.
"""
function zoom_interp(psi,x,y,xv,yv;win=1,nz=30)
    psiz, xz, yz = zoom_grid(psi,x,y,xv,yv,win=win)
    psi_itp = interpolate((xz,yz), psiz, Gridded(Linear()))
    xint,yint = LinRange(xz[1],xz[end],nz),LinRange(yz[1],yz[end],nz)
    psi_int = psi_itp(xint,yint)
    return psi_int,xint,yint
end

"""
    vortz,psiz = corezoom(vortex,ψ<:Field,winhalf=2,Nz=30)

Uses local interpolation to resolve core location.
"""
function corezoom(vortex::PointVortex,psi::T,winhalf=2,Nz=30) where T<:Field
    @unpack ψ,x,y = psi
    xv,yv,qv = vortex_array(vortex)
    dx,dy = Δ(x),Δ(y)
    ixv = isapprox.(x,xv,atol=dx) |> findfirst
    iyv = isapprox.(y,yv,atol=dy) |> findfirst
    ixwin = (ixv-winhalf):(ixv+winhalf-1)
    iywin = (iyv-winhalf):(iyv+winhalf-1)
    xw = x[ixwin]; yw = y[iywin]; psiw = ψ[ixwin,iywin]
    knots = (xw,yw)
    itp = interpolate(knots, psiw, Gridded(Linear()))
    xz = LinRange(xw[1],xw[end],Nz)
    yz = LinRange(yw[1],yw[end],Nz)
    psiz = itp(xz,yz)
    ψv = T(psiz,xz |> Vector,yz |> Vector)
    vortz = findvortices_grid(ψv,shift=true)
    vortz = remove_vortices_edge(vortz,ψv)[1]
    return vortz,ψv
end

corezoom(vortex::Array{PointVortex,1},psi::Field,winhalf=2,Nz=30) = corezoom(vortex[1],psi,winhalf,Nz)

function corezoom_periodic(vortex::PointVortex,psi::T,winhalf=2,Nz=30) where T<:Field
    @unpack ψ,x,y = psi
    xv,yv,qv = vortex_array(vortex)
    dx,dy = Δ(x),Δ(y)
    nx,ny = size(ψ)
    ixv = isapprox.(x,xv,atol=dx) |> findfirst
    iyv = isapprox.(y,yv,atol=dy) |> findfirst
    ix1,ix2,iy1,iy2 = ixv-winhalf,ixv+winhalf-1,iyv-winhalf,iyv+winhalf-1
    ixwin,iywin = ix1:ix2,iy1:iy2
    x1,x2 = x[ixv]-dx*winhalf,x[ixv]+dx*(winhalf-1)
    y1,y2 = y[iyv]-dy*winhalf,y[iyv]+dy*(winhalf-1)
    # periodic patch
    ixwinp = mod1.(ixwin,nx) # periodic indices
    iywinp = mod1.(iywin,ny)
    xw = x1:dx:x2   # window of x values
    yw = y1:dy:y2   # window of y values
    psiw = ψ[ixwinp,iywinp] # window of wavefunction
    itp = interpolate((xw,yw), psiw, Gridded(Linear()))
    xz,yz = LinRange(x1,x2,Nz),LinRange(y1,y2,Nz)
    ψv = T(itp(xz,yz),xz |> Vector,yz |> Vector)
    vortz = findvortices_grid(ψv,shift=true)
    vortz = remove_vortices_edge(vortz,ψv)[1]
    return vortz,ψv
end

corezoom_periodic(vortex::Array{PointVortex,1},psi::Field,winhalf=2,Nz=30) = corezoom_periodic(vortex[1],psi,winhalf,Nz)
