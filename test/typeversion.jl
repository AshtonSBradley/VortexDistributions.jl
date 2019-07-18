using Parameters, Interpolations, Plots, Test

abstract type FieldTopology end
abstract type Vortex end
abstract type VortexCore end

struct Ansatz <: VortexCore end
struct Exact <: VortexCore end

@with_kw mutable struct Torus <: FieldTopology
    ψ::Array{Complex{Float64},2}
    x::Vector{Float64}
    y::Vector{Float64}
end

@with_kw mutable struct Sphere <: FieldTopology
    ψ::Array{Complex{Float64},2}
    x::Vector{Float64}
    y::Vector{Float64}
end

@with_kw mutable struct PointVortex <: Vortex
    xv::Float64
    yv::Float64
    qv::Int64
end

PointVortex(v::Array{Float64,2}) = PointVortex.(v[:,1],v[:,2],v[:,3])
RawVortex(v::PointVortex) = [v.xv v.yv v.qv]
RawVortex(v::Array{PointVortex,1}) = reduce(vcat,RawVortex.(v))
randvortex(n) = PointVortex.(randn(n), randn(n), rand([1,-1],n))
function randvortex(n,psi::FieldTopology)
    @unpack ψ,x,y = psi
    return PointVortex.(rand(x,n),rand(y,n),rand([1,-1],n))
end

include("../src/detection.jl")

function findwhere(A)
    I = findall(!iszero,A)
    v = A[I]
    ix = [I[i][1] for i in eachindex(I)]
    iy = [I[i][2] for i in eachindex(I)]
    return ix,iy,v
end

function findvortices_jumps(psi::FieldTopology;shift=true)
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

    nt = np + nn

    if shift
        dx = x[2]-x[1]; dy = y[2] - y[1]
        xp .+= -dx/2; yp .+= -dy/2; xn .+= -dx/2; yn .+= -dy/2
    end
    vortices = [xn yn -vn; xp yp vp]

    return PointVortex(vortices)
end

Λ = 0.8249
r(x,y) = sqrt(x^2+y^2)
vcore(r) = (Λ*r)^2/(1+(Λ*r)^2)
vdensity(x,y) = sqrt(vcore(r(x,y)))
vphase(x,y) = atan(y,x)
vortexgrid(x,y,x0,y0,q0) = vdensity(x - x0, y - y0)*exp(im*q0*vphase(x-x0,y-y0))

function makevortex!(psi::Torus,vort::PointVortex)
    @unpack ψ,x,y = psi
    xv,yv,qv = RawVortex(vort)
    @. ψ *= vortexgrid(x,y',xv,yv,qv)
    @pack! psi = ψ
    return psi
end

function makevortex!(psi::Torus,vort::Array{PointVortex,1})
    for j in eachindex(vort)
        makevortex!(psi,vort[j])
    end
end

function findvortices_grid(psi::Torus;shift=true)
    vort = findvortices_jumps(psi,shift=shift)
    rawvort = RawVortex(vort)
    vort = sortslices(rawvort,dims=1)
    return PointVortex(vort)
end

function findvortices_grid(psi::Sphere;shift=true)
    @unpack ψ = psi
    windvals = phasejumps(angle.(ψ),2)
    vort = findvortices_jumps(psi,shift=shift)
    rawvort = RawVortex(vort)

    # North pole winding. - sign means polar vortex co-rotating with w > 0
    w1 = sum(windvals[1,:])
    (sign(w1) != 0) && (rawvort = [rawvort; 0.0 0.0 -w1])
    # South pole winding. + sign means co-rotating with w > 0
    w2 = sum(windvals[end,:])
    (sign(w2) != 0) && (rawvort = [rawvort; pi 0.0 w2])

    vort = sortslices(vort,dims=1)
    return PointVortex(vort)
end


function findvortices_interp(psi::FieldTopology)
    vort = findvortices_grid(psi,shift=false)
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
    vortices = findvortices(psi::T,interp=true) where T<:FieldTopology

Locates vortices as 2π phase windings around plaquettes on a cartesian spatial grid. Uses an optimized plaquette method followed by recursive interpolation.

Requires a 2D wavefunction ψ(x,y) on a cartesian grid specified by vectors x, y.

`vortices` - array of vortex coordinates `xv,yv` and charges `qv`. Each row is of the form `[xv, yv, cv]`, and the array is sorted into lexical order according to the `xv` coordinates
"""
function findvortices(psi::FieldTopology,interp=true)
    if interp
        return findvortices_interp(psi)
    else
        return findvortices_grid(psi)
    end
end

"""
`vortices = removeedgevortices(vort::Array{PointVortex,1},x,y,edge=1)`

Strips edgevortices due to periodic phase differencing."""
function removeedgevortices(vort::Array{PointVortex,1},psi::FieldTopology,edge=1)
    @unpack x,y = psi; dx=x[2]-x[1]; dy=y[2]-y[1]
    keep = []
    for j = 1:length(vort)
        xi,yi,qi = RawVortex(vort[j])
        xedge = isapprox(xi,x[1],atol=edge*dx) || isapprox(xi,x[end],atol=edge*dx)
        yedge = isapprox(yi,y[1],atol=edge*dy) || isapprox(yi,y[end],atol=edge*dy)
        not_edge = !(xedge || yedge)
        not_edge && push!(keep,j)
    end
    return vort[keep]
end

"""
    vortz,psiz,xz,yz = corezoom(vortex,psi::T,winhalf=2,Nz=30) where T<: FieldTopology

Uses local interpolation to resolve core location to ~ 5 figures.
"""
function corezoom(vortex::PointVortex,psi::T,winhalf=2,Nz=30) where T<:FieldTopology
    @unpack ψ,x,y = psi
    xv,yv,qv = RawVortex(vortex)
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
    vortz = findvortices_grid(ψv,shift=false) #TODO Really need this???
    vortz = removeedgevortices(vortz,ψv)
    return vortz[1],ψv
end

corezoom(vortex::Array{PointVortex,1},psi::FieldTopology,winhalf=2,Nz=30) = corezoom(vortex[1],psi,2,30)

N = 100
L = 100.0
x = LinRange(-L/2,L/2,N)
y = x
psi0 = one.(x*y') |> complex
psi = Torus(psi0,x,y)

vtest = randvortex(500,psi)
makevortex!(psi,vtest)
@unpack ψ = psi
heatmap(x,y,angle.(ψ),transpose=true)

vgridraw = findvortices_grid(psi,shift=false)
vgridraw = removeedgevortices(vgridraw,psi)
vgrid = findvortices_grid(psi)
vgrid = removeedgevortices(vgrid,psi)
vint = findvortices_interp(psi)

# ======================
# ==== Masking
function findvortexmask(psi::T,R) where T<: FieldTopology
    @unpack x,y,ψ = psi
    ψ = circmask(ψ,x,y,1.1*R)
    ϕ = T(ψ,x,y)
    vortices = findvortices(ϕ)

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
