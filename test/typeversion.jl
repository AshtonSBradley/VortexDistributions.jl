using Plots, JLD2, Parameters

@load "./examples/one_frame.jld2"

gr(transpose=false)
heatmap(abs2.(ψ1'))
heatmap(x,y,angle.(ψ1'))

xind = 100:400 #100:300 for no vortices
yind = 200:800
xwin = x[xind]
ywin = y[yind]
ψwin = ψ1[xind,yind]
heatmap(xwin,ywin,abs2.(ψwin'))
heatmap(xwin,ywin,angle.(ψwin'))

abstract type Topology end

struct Torus <: Topology
    x::Vector{Float64}
    y::Vector{Float64}
    ψ::Array{Complex{Float64},2}
end

struct Sphere <: Topology
    x::Vector{Float64}
    y::Vector{Float64}
    ψ::Array{Complex{Float64},2}
end

psi = Torus(x,y,ψ1)
psiw = Torus(xwin,ywin,ψwin)

function findallwhere(A)
    I = findall(!iszero,A)
    v = A[I]
    ix = [I[i][1] for i in eachindex(I)]
    iy = [I[i][2] for i in eachindex(I)]
    return ix,iy,v
end

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

function findvortices(psi::Torus)
    @unpack x,y,ψ = psi
    phase = angle.(ψ)
    diffx = phasejumps(phase,1)
    diffy = phasejumps(phase,2)

    #use in-place memory recycling
    circshift!(phase,diffx,(0,1))
    diffx .-= phase
    diffx .-= diffy
    circshift!(phase,diffy,(1,0))
    diffx .+= phase

    ixp,iyp,vp = findallwhere(diffx.>0.)
    xp = x[ixp]; yp = y[iyp]; np = length(vp)

    ixn,iyn,vn = findallwhere(diffx.<0.)
    xn = x[ixn]; yn = y[iyn]; nn = length(vn)

    nt = np + nn

    dx = x[2]-x[1]; dy = y[2]-y[1]
    xp .+= -dx/2; yp .+= -dy/2; xn .+= -dx/2; yn .+= -dy/2

    vortices = [xn yn -vn; xp yp vp]
    vortices = sortslices(vortices,dims=1)

    return nt,np,nn,vortices
end

function findvortices(psi::Sphere)
    @unpack x,y,ψ = psi
    phase = angle.(ψ)
    diffx = phasejumps(phase,1)
    diffy = phasejumps(phase,2)

    #use in-place memory recycling
    circshift!(phase,diffx,(0,1))
    diffx .-= phase
    diffx .-= diffy
    circshift!(phase,diffy,(1,0))
    diffx .+= phase

    ixp,iyp,vp = findallwhere(diffx.>0.)
    xp = x[ixp]; yp = y[iyp]; np = length(vp)

    ixn,iyn,vn = findallwhere(diffx.<0.)
    xn = x[ixn]; yn = y[iyn]; nn = length(vn)

    nt = np + nn

    dx = x[2]-x[1]; dy = y[2] - y[1]
    xp .+= -dx/2; yp .+= -dy/2; xn .+= -dx/2; yn .+= -dy/2

    #find polar winding
    vortices = [xn yn -vn; xp yp vp]

    # to do: optimize this step at start of method to avoid extra angle call:
    windvals = phasejumps(angle.(ψ),2)

    #test north pole for winding. - sign means polar vortex co-rotating with w > 0
    w1 = sum(windvals[1,:])
    (sign(w1) != 0) && (vortices = [vortices; 0.0 0.0 -w1])

    #test south pole for winding. + sign means co-rotating with w > 0
    w2 = sum(windvals[end,:])
    (sign(w2) != 0) && (vortices = [vortices; pi 0.0 w2])

    vortices = sortslices(vortices,dims=1)
    return nt,np,nn,vortices
end

#this test appears to work:
nt,np,nn,vortices = findvortices(psiw)


function findvortices_interp(psi,x,y)
    nt,np,nn,vortices = findvortices_grid(psi,x,y)
    nt,np,nn,vortices = remove_edgevortices(vortices,x,y)
    #TODO: allow for interp with periodic data (here edges are stripped to avoid)
    for j in 1:nt
        vortex = vortices[j,:]
        vortz,psiz,xz,yz = corezoom(vortex,psi,x,y)
        vortz,psiz,xz,yz = corezoom(vortz,psiz,xz,yz)
        vortz,psiz,xz,yz = corezoom(vortz,psiz,xz,yz)
        vortices[j,1:2] = vortz[1:2]
    end
    return nt,np,nn,vortices
end

"""
   `nt,np,nn,vortices = findvortices(ψ,x,y)`

Locates vortices as 2π phase windings around plaquettes on a cartesian spatial grid. Uses an optimized plaquette method followed by recursive interpolation.

Requires a 2D wavefunction ψ(x,y) on a cartesian grid specified by vectors x, y.

`nt` - total number of vortices

`np` - number of positive vortices

`nn` - number of negative vortices

`vortices` - array of vortex coordinates `xv,yv` and circulations `cv`. Each row is of the form `[xv, yv, cv]`, and the array is sorted into lexical order according to the `xv` coordinates
"""
function findvortices(psi,x,y,interp::Bool=true)
    if interp
        nt,np,nn,vortices = findvortices_interp(psi,x,y)
    elseif !interp
        nt,np,nn,vortices = findvortices_grid(psi,x,y)
    end
    return nt,np,nn,vortices
end
