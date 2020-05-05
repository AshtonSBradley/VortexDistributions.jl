
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
    end
    return vort
end

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

"""
    vortz,psiz = corezoom(vortex,ψ<:Field,winhalf=2,Nz=30)

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


## phase jumps, unwrap

"""
    jumps = phasejumps(ϕ,dim)

Count phase jumps greater than `π` in phase `ϕ` along dimension `dim`.
See also: [`unwrap`](@ref), [`unwrap`](@ref), [`unwrap!`](@ref)
"""
function phasejumps(phase,dim=1)
    @assert (dim==1 || dim==2)
    s1,s2 = size(phase)
    pdiff = zero(phase)

    if dim == 1
    for j in 1:s2
        @inbounds Δ = phase[1,j] - phase[s1,j]
        @inbounds abs(Δ) > π && (pdiff[1,j] += sign(Δ))
        for i in 2:s1
            @inbounds Δ = phase[i,j] - phase[i-1,j]
            @inbounds abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end

    end

    elseif dim == 2
    for i in 1:s1
        @inbounds Δ = phase[i,1] - phase[i,s2]
        @inbounds abs(Δ) > π && (pdiff[i,1] += sign(Δ))
        for j in 2:s2
            @inbounds Δ = phase[i,j] - phase[i,j-1]
            @inbounds abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end
    end
    end
  return pdiff
end

function phasejumps!(pdiff,phase,dim=1)
    @assert (dim==1 || dim==2)
    s1,s2 = size(phase)

    if dim == 1
    for j in 1:s2
        @inbounds Δ = phase[1,j] - phase[s1,j]
        @inbounds abs(Δ) > π && (pdiff[1,j] += sign(Δ))
        for i in 2:s1
            @inbounds Δ = phase[i,j] - phase[i-1,j]
            @inbounds abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end

    end

    elseif dim == 2
    for i in 1:s1
        @inbounds Δ = phase[i,1] - phase[i,s2]
        @inbounds abs(Δ) > π && (pdiff[i,1] += sign(Δ))
        for j in 2:s2
            @inbounds Δ = phase[i,j] - phase[i,j-1]
            @inbounds abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end
    end
    end
end

"""
    ϕu = unwrap(ϕ,dim=1)

Unwraps 2d array `ϕ` along dimension `dim`, acting periodically to give back array of same size as `ϕ`.
See also: [`unwrap!`](@ref)
"""
function unwrap(phase::Array{Float64,2},dim=1)
    @assert (dim==1 || dim==2)
    s1,s2 = size(phase)
    uphase = copy(phase)

    if dim == 1
    for j in 1:s2
        @inbounds uphase[1,j] - uphase[s1,j] >= π && (uphase[1,j] -= 2π)
        @inbounds uphase[1,j] - uphase[s1,j] <= -π && (uphase[1,j] += 2π)
        for i in 2:s1
        @inbounds uphase[i,j] - uphase[i-1,j] >= π && (uphase[i,j] -= 2π)
        @inbounds uphase[i,j] - uphase[i-1,j] <= -π && (uphase[i,j] += 2π)
        end
    end

    elseif dim == 2
    for i in 1:s1
            @inbounds uphase[i,1] - uphase[i,s2] >= π && (uphase[i,1] -= 2π)
            @inbounds uphase[i,1] - uphase[i,s2] <= -π && (uphase[i,1] += 2π)
        for j in 2:s2
            @inbounds uphase[i,j] - uphase[i,j-1] >= π && (uphase[i,j] -= 2π)
            @inbounds uphase[i,j] - uphase[i,j-1] <= -π && (uphase[i,j] += 2π)
        end
    end
end

    return uphase
end

"""
    unwrap!(ϕu,ϕ,dim)

Unwraps 2d phase array `ϕ` along dimension `dim`,
acting periodically and writing the unwrapped array `ϕ`
in place.
See also: [`unwrap`](@ref)
"""
function unwrap!(uphase::Array{Float64,2},phase::Array{Float64,2},dim=1)
    @assert (dim==1 || dim==2)
    s1,s2 = size(phase)
    uphase .= phase

    if dim == 1
    for j in 1:s2
        @inbounds uphase[1,j] - uphase[s1,j] >= π && (uphase[1,j] -= 2π)
        @inbounds uphase[1,j] - uphase[s1,j] <= -π && (uphase[1,j] += 2π)
        for i in 2:s1
        @inbounds uphase[i,j] - uphase[i-1,j] >= π && (uphase[i,j] -= 2π)
        @inbounds uphase[i,j] - uphase[i-1,j] <= -π && (uphase[i,j] += 2π)
        end
    end

    elseif dim == 2
    for i in 1:s1
            @inbounds uphase[i,1] - uphase[i,s2] >= π && (uphase[i,1] -= 2π)
            @inbounds uphase[i,1] - uphase[i,s2] <= -π && (uphase[i,1] += 2π)
        for j in 2:s2
            @inbounds uphase[i,j] - uphase[i,j-1] >= π && (uphase[i,j] -= 2π)
            @inbounds uphase[i,j] - uphase[i,j-1] <= -π && (uphase[i,j] += 2π)
        end
    end
end
end

function unwrap(phase::Array{Float64,1})
    uphase = copy(phase)
    s1 = length(phase)
        @inbounds uphase[1] - uphase[s1] >= π && (uphase[1] -= 2π)
        @inbounds uphase[1] - uphase[s1] <= -π && (uphase[1] += 2π)
        for i in 2:s1
        @inbounds uphase[i] - uphase[i-1] >= π && (uphase[i] -= 2π)
        @inbounds uphase[i] - uphase[i-1] <= -π && (uphase[i] += 2π)
        end

        return uphase
end

function unwrap!(uphase::Array{Float64,1},phase::Array{Float64,1})
    s1 = length(phase)
        @inbounds uphase[1] - uphase[s1] >= π && (uphase[1] -= 2π)
        @inbounds uphase[1] - uphase[s1] <= -π && (uphase[1] += 2π)
        for i in 2:s1
        @inbounds uphase[i] - uphase[i-1] >= π && (uphase[i] -= 2π)
        @inbounds uphase[i] - uphase[i-1] <= -π && (uphase[i] += 2π)
        end
end
