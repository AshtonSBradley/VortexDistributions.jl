using Parameters, Interpolations

abstract type FieldTopology end
abstract type Vortex end
abstract type VortexCore end

struct Ansatz <: VortexCore end
struct Exact <: VortexCore end

struct Torus <: FieldTopology
    x::Vector{Float64}
    y::Vector{Float64}
    ψ::Array{Complex{Float64},2}
end

struct Sphere <: FieldTopology
    x::Vector{Float64}
    y::Vector{Float64}
    ψ::Array{Complex{Float64},2}
end

struct PointVortex <: Vortex
    xv::Float64
    yv::Float64
    qv::Int64
end

struct Vortex3 <: Vortex
    x::Float64
    y::Float64
    z::Float64
    q::Int64
end

function find_on_grid_unsorted(psi::FieldTopology,shift=true) where T<:FieldTopology
    @unpack x,y,ψ = psi
    phase = angle.(ψ)
    diffx = phasejumps(phase,1); diffy = phasejumps(phase,2)

    # use in-place memory recycling
    circshift!(phase,diffx,(0,1))
    diffx .-= phase; diffx .-= diffy
    circshift!(phase,diffy,(1,0))
    diffx .+= phase

    ixp,iyp,vp = findwhere(diffx .> 0.0)
    xp = x[ixp]; yp = y[iyp]; np = length(vp)

    ixn,iyn,vn = findwhere(diffx .< 0.0)
    xn = x[ixn]; yn = y[iyn]; nn = length(vn)

    nt = np + nn

    if shift
        dx = x[2]-x[1]; dy = y[2] - y[1]
        xp .+= -dx/2; yp .+= -dy/2; xn .+= -dx/2; yn .+= -dy/2
    end
    vortices = [xn yn -vn; xp yp vp]

    return nt,np,nn,vortices
end

function findvorticesgrid(psi::Torus,shift=true)
    @unpack x,y,ψ = psi

    nt,np,nn,vortices = find_on_grid_unsorted(psi,shift)

    vortices = sortslices(vortices,dims=1)

    return nt,np,nn,vortices
end

function findvorticesgrid(psi::Sphere,shift=true)

    nt,np,nn,vortices = find_on_grid_unsorted(psi,shift)

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


function findvorticesinterp(psi::T) where T<:FieldTopology
    @unpack ψ,x,y = psi
    nt,np,nn,vortices = findvorticesgrid(psi;shift=false)
    nt,np,nn,vortices = removeedgevortices(vortices,x,y)
    #TODO: allow for interp with periodic data (here edges are stripped to avoid)
    for j in 1:nt
        vortex = vortices[j,:]
        vortz,psiz,xz,yz = corezoom(vortex,ψ,x,y)
        vortz,psiz,xz,yz = corezoom(vortz,psiz,xz,yz)
        vortz,psiz,xz,yz = corezoom(vortz,psiz,xz,yz)
        vortices[j,1:2] = vortz[1:2]
    end
    return nt,np,nn,vortices
end

"""
    nt,np,nn,vortices = findvortices(psi::T,interp=true) where T<:FieldTopology

Locates vortices as 2π phase windings around plaquettes on a cartesian spatial grid. Uses an optimized plaquette method followed by recursive interpolation.

Requires a 2D wavefunction ψ(x,y) on a cartesian grid specified by vectors x, y.

`nt` - total number of vortices

`np` - number of positive vortices

`nn` - number of negative vortices

`vortices` - array of vortex coordinates `xv,yv` and circulations `cv`. Each row is of the form `[xv, yv, cv]`, and the array is sorted into lexical order according to the `xv` coordinates
"""
function findvortices(psi::T;interp=true) where T<:FieldTopology
    if interp
        return nt,np,nn,vortices = findvorticesinterp(psi)
    else
        return nt,np,nn,vortices = findvorticesgrid(psi)
    end
end

function findvortexmask(psi::T,R) where T<: FieldTopology
    @unpack x,y,ψ = psi
    ψ = circmask(ψ,x,y,1.1*R)
    ϕ = T(x,y,ψ)
    nt,np,nn,vortices = findvortices(ϕ)

# TODO more tests
# remove vortices found outside mask boundary
for (i,xv) in enumerate(vortices[:,1]), yv in vortices[i,2]
    (norm([xv,yv]) > R) && (vortices[i,:] = [0. 0. 0.])
end

vortfound = Array{Float64, 2}[]
for (i,sv) in enumerate(vortices[:,3])
    (sv!=0) && push!(vortfound,vortices[i,:]')
end
    # currently will fail for more that one vortex (?)
return vortfound[1]
end

"""
`vortices = remove_edgevortices(vortices,x,y,edge=1)`

Strips edgevortices due to periodic phase differencing."""
function removeedgevortices(vortices,x,y,edge=1)
    dx=x[2]-x[1]; dy=y[2]-y[1]
    Nfound,_ = size(vortices)
    keep = []
    for j = 1:Nfound
        xi = vortices[j,1]; yi = vortices[j,2];
        xedge = isapprox(xi,x[1],atol=edge*dx) || isapprox(xi,x[end],atol=edge*dx)
        yedge = isapprox(yi,y[1],atol=edge*dy) || isapprox(yi,y[end],atol=edge*dy)
        not_edge = !(xedge || yedge)
        not_edge && push!(keep,j)
    end
    vortices = vortices[keep,:]
    nt = size(vortices)[1]
    np = sum(vortices[:,3].>0)
    nn = nt - np
    return nt,np,nn,vortices
end

"""
    vortz,psiz,xz,yz = corezoom(vortices,psi::T,winhalf=2,Nz=30) where T<: FieldTopology

Uses local interpolation to resolve core location to ~ 5 figures.
"""
function corezoom(vortices,psi::T,winhalf=2,Nz=30) where T<:FieldTopology
    @unpack ψ,x,y = psi
    xv,yv = vortices[1:2]
    dx=x[2]-x[1];dy=y[2]-y[1]
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
    ψv = T(xz |> Vector,yz |> Vector,psiz)
    nt,np,nn,vortz = findvorticesgrid(ψv;shift=false)
    nt,np,nn,vortz = removeedgevortices(vortz,xz |> Vector,yz |> Vector)
    return vortz,psiz,xz,yz
end


#some helpers
linspace(a,b,n) = LinRange(a,b,n) |> collect

function edgemask(psi::T) where T<:FieldTopology
    @unpack ψ,x,y = psi
    ψ[:,1] = zero(x)
    ψ[:,end] = zero(x)
    ψ[1,:] = zero(y')
    ψ[end,:] = zero(y')
    return T(x,y,ψ)
end

function circmask(psi::T,R) where T<:FieldTopology
    @unpack ψ,x,y = psi
    for j in eachindex(x), k in eachindex(y)
        (x[j]^2+y[k]^2 > R^2) && (ψ[j,k] = complex(0.))
    end
    return ψ
end

function circmask!(phi::T,psi::T,R) where T<: FieldTopology
    @unpack ψ,x,y = psi
    phi.ψ .= psi.ψ
    for j in eachindex(x), k in eachindex(y)
            (x[j]^2+y[k]^2 > R^2) && (phi.ψ[j,k] = complex(0.))
    end
    return nothing
end

function isinterior(a,b,x,y)
    Lx = x[end]-x[1]; Ly = y[end] - y[1]
    dx = x[2] - x[1]; dy = y[2] - y[1]
    return (-Lx/2 + 2dx < a < Lx/2 - 2dx && -Ly/2 + 2dy < b < Ly/2 - 2dy)
end

function randomvortices(x,y,Nv)
    Lx = x[end]-x[1]; Ly = y[end] - y[1]
    dx = x[2] - x[1]; dy = y[2] - y[1]
    testvort = zeros(Nv,3)

    k = 1
    while k<=Nv
        a = -Lx/2 + Lx*rand()
        b = -Ly/2 + Ly*rand()
        σ = rand([-1 1])

        #make sure vortices are away from edges
        if isinterior(a,b,x,y)
            testvort[k,:] = [a b σ]
            k+=1
        end
    end
    return sortslices(testvort,dims=1)
end

function checkvortexlocations(testvort,vortices,x,y,Nv)
#check detection to 2 x grid resolution
    dx = x[2] - x[1]; dy = y[2] - y[1]
    vortfound = 0
    for j=1:Nv
        if (isapprox(testvort[j,1],vortices[j,1],atol=2dx) && isapprox(testvort[j,2],vortices[j,2],atol=2dy))
        vortfound+=1
        end
    end
    return vortfound
end



function makevortex(ψ,vortex,x,y,ξ=1.0)
    @assert typeof(x)==Array{Float64,1}
    @assert typeof(y)==Array{Float64,1}
    @assert typeof(ψ)==Array{Complex{Float64},2}
    x0, y0, σ0 = vortex
    R(x,y) = sqrt(x^2+y^2)
    return  @. ψ*vortexcore(R(x.-x0,y'.-y0),ξ)*exp(im*σ0*atan(y'.-y0,x.-x0))
end

function makevortex!(ψ,vortex,x,y,ξ=1.0)
    @assert typeof(x)==Array{Float64,1}
    @assert typeof(y)==Array{Float64,1}
    @assert typeof(ψ)==Array{Complex{Float64},2}
    x0, y0, σ0 = vortex
    R(x,y) = sqrt(x^2+y^2)
    ψ .= @. ψ*vortexcore(R(x.-x0,y'.-y0),ξ)*exp(im*σ0*atan(y'.-y0,x.-x0))
end

function makeallvortices!(ψ,vortices,x,y,ξ=1.0)
    for j = 1:size(vortices)[1]
        makevortex!(ψ,vortices[j,:],x,y,ξ)
    end
end



vort = PointVortex(.1,.2,1)

x = randn(10);y = randn(10); charge = rand([-1,1])

vort = PointVortex.(x,y,charge)

vort[3]
