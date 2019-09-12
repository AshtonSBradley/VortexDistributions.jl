function foundNear(n)
    near = true
    for j ∈ 1:n
        psi,vort = randVortexField(1)
        vortfound = findvortices(psi)
        vfdata = rawData(vortfound)
        vdata = rawData(vort)
        near *= isapprox(vdata,vfdata,rtol = 0.2)
    end
    return near
end
function findwhere(A)
    I = findall(!iszero,A)
    v = A[I]
    ix = [I[i][1] for i in eachindex(I)]
    iy = [I[i][2] for i in eachindex(I)]
    return ix,iy,v
end

function findvortices_jumps(psi::Field;shift=true)
    @unpack x,y,ψ = psi
    phase = angle.(ψ)

    # Compute all nearest neighbour jumps
    diffx = phasejumps(phase,1); diffy = phasejumps(phase,2)

    # Compute all plaquette loop integrals with in-place memory recycling
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
        dx = x[2]-x[1]; dy = y[2] - y[1]
        xp .-= dx/2; yp .-= dy/2; xn .-= dx/2; yn .-= dy/2
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
    vort = removeedgevortices(vort,psi)
    #TODO: allow for interp with periodic data (here edges are stripped)

    for (j,vortex) in enumerate(vort)
        try
        vortz,psiz = corezoom(vortex,psi)
        vortz,psiz = corezoom(vortz,psiz)
        vortz,psiz = corezoom(vortz,psiz)
        vort[j] = vortz
        catch nothing
        end
        return vort
    end
end

"""
    vortices = findvortices(psi::T,interp=true) where T<:Field

Locates vortices as 2π phase windings around plaquettes on a cartesian spatial grid. Uses an optimized plaquette method followed by recursive interpolation.

Requires a 2D wavefunction ψ(x,y) on a cartesian grid specified by vectors x, y.

`vortices` - array of vortex coordinates `xv,yv` and charges `qv`. Each row is of the form `[xv, yv, cv]`, and the array is sorted into lexical order according to the `xv` coordinates
"""
function findvortices(psi::Field,interp=true)
    if interp
        return findvortices_interp(psi)
    else
        return findvortices_grid(psi)
    end
end


"""
`vortices = removeedgevortices(vort::Array{PointVortex,1},x,y,edge=1)`

Strips edgevortices due to periodic phase differencing."""
function removeedgevortices(vort::Array{PointVortex,1},psi::Field,edge=1)
    @unpack x,y = psi; dx=x[2]-x[1]; dy=y[2]-y[1]
    keep = []
    for j = 1:length(vort)
        xi,yi,qi = rawData(vort[j])
        xedge = isapprox(xi,x[1],atol=edge*dx) || isapprox(xi,x[end],atol=edge*dx)
        yedge = isapprox(yi,y[1],atol=edge*dy) || isapprox(yi,y[end],atol=edge*dy)
        not_edge = !(xedge || yedge)
        not_edge && push!(keep,j)
    end
    return vort[keep]
end

"""
    vortz,psiz = corezoom(vortex,psi::T,winhalf=2,Nz=30) where T <: Field

Uses local interpolation to resolve core location.
"""
function corezoom(vortex::PointVortex,psi::T,winhalf=2,Nz=30) where T<:Field
    @unpack ψ,x,y = psi
    xv,yv,qv = rawData(vortex)
    dx=x[2]-x[1]; dy=y[2]-y[1]
    ixv = isapprox.(x,xv,atol=dx) |> findlast
    iyv = isapprox.(y,yv,atol=dy) |> findlast
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


"""
    jumps = phasejumps(phase,dim)

Count phase jumps greater than π in `phase` along dimension `dim`
"""
function phasejumps(phase,dim=1)
    @assert (dim==1 || dim==2)
    Nx,Ny = size(phase)
    pdiff = zero(phase)

    if dim == 1
    @inbounds for j in 1:Ny
        @inbounds for i in 2:Nx
            Δ = phase[i,j] - phase[i-1,j]
            abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end
            Δ = phase[1,j] - phase[Nx,j]
            abs(Δ) > π && (pdiff[1,j] += sign(Δ))
    end

    elseif dim == 2
    @inbounds for j in 2:Ny
        @inbounds for i in 1:Nx
            Δ = phase[i,j] - phase[i,j-1]
            abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end
    end
        @inbounds for i in 1:Nx
            Δ = phase[i,1] - phase[i,Ny]
            abs(Δ) > π && (pdiff[i,1] += sign(Δ))
        end
    end
  return pdiff
end

function phasejumps!(pdiff,phase,dim=1)
    @assert (dim==1 || dim==2)
    Nx,Ny = size(phase)

    if dim == 1
    @inbounds for j in 1:Ny
        @inbounds for i in 2:Nx
            Δ = phase[i,j] - phase[i-1,j]
            abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end
            Δ = phase[1,j] - phase[Nx,j]
            abs(Δ) > π && (pdiff[1,j] += sign(Δ))
    end

    elseif dim == 2
    @inbounds for j in 2:Ny
        @inbounds for i in 1:Nx
            Δ = phase[i,j] - phase[i,j-1]
            abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end

    end
        @inbounds for i in 1:Nx
            Δ = phase[i,1] - phase[i,Ny]
            abs(Δ) > π && (pdiff[i,1] += sign(Δ))
        end
    end
end

"""
`unwrapped = unwrap(phase,dim=1)`

Unwraps 2d array `phase` along dimension `dim`, acting periodically to give back array of same size as `phase`.

`unwrap!(unwrapped,phase,dim)` writes in-place to `unwrapped`.
"""
function unwrap(phase::Array{Float64,2},dim=1)
    @assert (dim==1 || dim==2)
    Nx,Ny = size(phase)
    uphase = copy(phase)

    if dim == 1
    @inbounds for j in 1:Ny
        @inbounds for i in 2:Nx
        (uphase[i,j] - uphase[i-1,j] >= π) && (uphase[i,j] -= 2π)
        (uphase[i,j] - uphase[i-1,j] <= -π) && (uphase[i,j] += 2π)
        end
        (uphase[1,j] - uphase[Nx,j] >= π) && (uphase[1,j] -= 2π)
        (uphase[1,j] - uphase[Nx,j] <= -π) && (uphase[1,j] += 2π)
    end

    elseif dim == 2
    @inbounds for j in 2:Ny
        @inbounds for i in 1:Nx
        (uphase[i,j] - uphase[i,j-1] >= π) && (uphase[i,j] -= 2π)
        (uphase[i,j] - uphase[i,j-1] <= -π) && (uphase[i,j] += 2π)
        end

    end
        @inbounds for i in 1:Nx
        (uphase[i,1] - uphase[i,Ny] >= π) && (uphase[i,1] -= 2π)
        (uphase[i,1] - uphase[i,Ny] <= -π) && (uphase[i,1] += 2π)
        end
    end

    return uphase
end

function unwrap!(uphase::Array{Float64,2},phase::Array{Float64,2},dim=1)
    @assert (dim==1 || dim==2)
    Nx,Ny = size(phase)
    uphase .= phase

    if dim == 1
    @inbounds for j in 1:Ny
        @inbounds for i in 2:Nx
        (uphase[i,j] - uphase[i-1,j] >= π) && (uphase[i,j] -= 2π)
        (uphase[i,j] - uphase[i-1,j] <= -π) && (uphase[i,j] += 2π)
        end
        (uphase[1,j] - uphase[Nx,j] >= π) && (uphase[1,j] -= 2π)
        (uphase[1,j] - uphase[Nx,j] <= -π) && (uphase[1,j] += 2π)
    end

    elseif dim == 2
    @inbounds for j in 2:Ny
        @inbounds for i in 1:Nx
        (uphase[i,j] - uphase[i,j-1] >= π) && (uphase[i,j] -= 2π)
        (uphase[i,j] - uphase[i,j-1] <= -π) && (uphase[i,j] += 2π)
        end

    end
        @inbounds for i in 1:Nx
        (uphase[i,1] - uphase[i,Ny] >= π) && (uphase[i,1] -= 2π)
        (uphase[i,1] - uphase[i,Ny] <= -π) && (uphase[i,1] += 2π)
        end
    end
end

function unwrap(phase::Array{Float64,1})
    uphase = copy(phase)
    Nx = length(phase)
        @inbounds for i in 2:Nx
        (uphase[i] - uphase[i-1] >= π) && (uphase[i] -= 2π)
        (uphase[i] - uphase[i-1] <= -π) && (uphase[i] += 2π)
        end
        (uphase[1] - uphase[Nx] >= π) && (uphase[1] -= 2π)
        (uphase[1] - uphase[Nx] <= -π) && (uphase[1] += 2π)
        return uphase
end

function unwrap!(uphase::Array{Float64,1},phase::Array{Float64,1})
    Nx = length(phase)
        @inbounds for i in 2:Nx
        (uphase[i] - uphase[i-1] >= π) && (uphase[i] -= 2π)
        (uphase[i] - uphase[i-1] <= -π) && (uphase[i] += 2π)
        end
        (uphase[1] - uphase[Nx] >= π) && (uphase[1] -= 2π)
        (uphase[1] - uphase[Nx] <= -π) && (uphase[1] += 2π)
end
