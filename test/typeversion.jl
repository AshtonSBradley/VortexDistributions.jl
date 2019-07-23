using Parameters, Interpolations, Plots, Test

abstract type Field end
abstract type Vortex end

@with_kw mutable struct Torus <: Field
    ψ::Array{Complex{Float64},2}
    x::Vector{Float64}
    y::Vector{Float64}
end

@with_kw mutable struct Sphere <: Field
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
RawData(v::PointVortex) = [v.xv v.yv v.qv]
RawData(v::Array{PointVortex,1}) = reduce(vcat,RawData.(v))

function uniform(a,b)
    @assert a<b
    return a + (b-a)*rand()
end

function uniform(a,b,n)
    a = a*ones(n)
    return uniform.(a,b)
end

uniform(n) = uniform(-.5,.5,n)
randcharge(n) = rand([-1 1],n)
randPointVortex(n) = PointVortex.(uniform(n), uniform(n), randcharge(n))

function randPointVortex(n,psi::Field)
    @unpack ψ,x,y = psi
    dx = x[2]-x[1]; dy = y[2]-y[1]
    xi,xf = x[1]+2dx,x[end]-2dx
    yi,yf = y[1]+2dy,y[end]-2dy
    return PointVortex.(uniform(xi,xf,n),uniform(yi,yf,n),randcharge(n))
end

#TODO
# function randvortex!(psi:Field,n)
#     @unpack ψ,x,y = psi
#     vort = randPointVortex(n,psi)
#     makevortex!()
# end

include("../src/detection.jl")

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

# Construction

# ansatz vortex wavefunction
Λ = 0.8249
r(x,y) = sqrt(x^2+y^2)
vcore(r) = sqrt((Λ*r)^2/(1+(Λ*r)^2))
vcore(x,y) = vcore(r(x,y))
vphase(x,y) = atan(y,x)
vortexgrid(x,y,x0,y0,q0) = vcore(x - x0, y - y0)*exp(im*q0*vphase(x-x0,y-y0))

function gpeansatz!(psi::Torus,vort::PointVortex)
    @unpack ψ,x,y = psi
    xv,yv,qv = RawData(vort)
    @. ψ *= vortexgrid(x,y',xv,yv,qv)
    @pack! psi = ψ
    return psi
end

function gpeansatz!(psi::Torus,vort::Array{PointVortex,1})
    for j in eachindex(vort)
        gpeansatz!(psi,vort[j])
    end
end

# TODO create and save some decent interpolations
# function gpeexact!(psi:Torus,vort::PointVortex)
#     @unpack ψ,x,y = psi
#     xv,yv,qv = RawData(vort)
#     if !isdefined(VortexDistributions,:gpecore)
#     loadpath = dirname(pathof(VortexDistributions))
#     @load joinpath(loadpath,"gpecore.jld2") y gpecore
#     end
#
#     return psi
# end

# function make_fastcore(n)
#     if n==1
#         loadpath = dirname(pathof(VortexDistributions))
#         @load loadpath*"/vortexcore.jld2" y ψ
#     else
#         y,ψ,res = gpecore(n)
#     end
#         ψi = interpolate(tuple(y[1:end-1]), ψ[1:end-1], Gridded(Linear()))
#     return ψi
# end

# Detection

function findvortices_grid(psi::Torus;shift=true)
    vort = findvortices_jumps(psi,shift=shift)
    rawvort = RawData(vort)
    vort = sortslices(rawvort,dims=1)
    return PointVortex(vort)
end

function findvortices_grid(psi::Sphere;shift=true)
    @unpack ψ = psi
    windvals = phasejumps(angle.(ψ),2)
    vort = findvortices_jumps(psi,shift=shift)
    rawvort = RawData(vort)

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
        xi,yi,qi = RawData(vort[j])
        xedge = isapprox(xi,x[1],atol=edge*dx) || isapprox(xi,x[end],atol=edge*dx)
        yedge = isapprox(yi,y[1],atol=edge*dy) || isapprox(yi,y[end],atol=edge*dy)
        not_edge = !(xedge || yedge)
        not_edge && push!(keep,j)
    end
    return vort[keep]
end

"""
    vortz,psiz,xz,yz = corezoom(vortex,psi::T,winhalf=2,Nz=30) where T<: Field

Uses local interpolation to resolve core location to ~ 5 figures.
"""
function corezoom(vortex::PointVortex,psi::T,winhalf=2,Nz=30) where T<:Field
    @unpack ψ,x,y = psi
    xv,yv,qv = RawData(vortex)
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

#TODO tidy this up!
corezoom(vortex::Array{PointVortex,1},psi::Field,winhalf=2,Nz=30) = corezoom(vortex[1],psi,2,30)

# Masking

function findvortexmask(psi::T,R) where T<: Field
    @unpack x,y,ψ = psi
    ψ = circmask(psi,1.1*R)
    ϕ = T(x,y,ψ)
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

# Helpers
linspace(a,b,n) = LinRange(a,b,n) |> collect

function edgemask!(psi::T) where T<:Field
    @unpack ψ,x,y = psi
    ψ[:,1] = zero(x) |> complex
    ψ[:,end] = zero(x) |> complex
    ψ[1,:] = zero(y') |> complex
    ψ[end,:] = zero(y') |> complex
    return T(ψ,x,y)
end

function circmask(psi::T,R) where T<:Field
    @unpack ψ,x,y = psi
    for j in eachindex(x), k in eachindex(y)
        (x[j]^2+y[k]^2 > R^2) && (ψ[j,k] = complex(0.))
    end
    return T(x,y,ψ)
end

function circmask!(phi::T,psi::T,R) where T<: Field
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

#TODO update this function
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


# Testing

N = 200
L = 100.0
x = LinRange(-L/2,L/2,N)
y = x
psi0 = one.(x*y') |> complex
psi = Torus(psi0,x,y)

vtest = randPointVortex(1,psi)
gpeansatz!(psi,vtest)
@unpack ψ = psi
heatmap(x,y,angle.(ψ),transpose=true)




# test detection
vgridraw = findvortices_grid(psi,shift=false)
vgridraw = removeedgevortices(vgridraw,psi)



vgrid = findvortices_grid(psi)
vgrid = removeedgevortices(vgrid,psi)



# creation
N = 200
L = 100.0
x = LinRange(-L/2,L/2,N)
y = x
psi0 = one.(x*y') |> complex
psi = Torus(psi0,x,y)

vtest = randPointVortex(1,psi)
xv,yv,qv = RawData(vtest)
y,ψ,res = gpecore(n)

function gpeexact!(psi:Torus,vort::PointVortex)
    @unpack ψ,x,y = psi
    xv,yv,qv = RawData(vort)
    if !isdefined(VortexDistributions,:gpecore)
    loadpath = dirname(pathof(VortexDistributions))
    @load joinpath(loadpath,"gpecore.jld2") y gpecore
    end

    return psi
end

"""
Make and evaluate the vortex core interpolation for charge `n`
"""
function core_chargen(x,n,ξ=1)
    y,ψ,res = gpecore(n)
    ψi = interpolate(tuple(y[1:end-1]), ψ[1:end-1], Gridded(Linear()))
    # isdefined(VortexDistributions,:ψi) && (ψi = make_fastcore(n))
    return ψi(x/ξ)
end

"""
Make the vortex core interpolation for charge `n`
"""
function make_fastcore(n)
    if n==1
        loadpath = dirname(pathof(VortexDistributions))
        @load loadpath*"/vortexcore.jld2" y ψ
    else
        y,ψ,res = gpecore(n)
    end
        ψi = interpolate(tuple(y[1:end-1]), ψ[1:end-1], Gridded(Linear()))
    return ψi
end

"""
Evaluate fast core interpolation
"""
function vortexcore(r,ψi::Interpolations.GriddedInterpolation,ξ=1)
    return ψi(r/ξ)
end

function gpecore(K,L=2,N=100,R = K)
    #currently r does  nothing!
    #N = 100
    #L = 2
    #Κ = 1
    #R = 2 # Stretching coordinate.  R ~ kappa seems to be a good ballpark

# Convert Chebyshev grid to [0,L] from [1 -1];
z,Dz = getChebDMatrix(N)
blank,D2z = getChebD2Matrix(N)
z = vec(reverse(z,dims=2))
z = (z)*L./2
Dz = -2*Dz./L
D2z = 4*D2z./L.^2

# y is the physical coordinate, range [0 infinity] (i.e y = r)
y =  R*(1 .+z)./(1 .-z)

#initial guess based on ansatz
Λ = 0.8248
ψ = @. y/sqrt(y^2 + Λ^-2)
ψ[1] = 0
ψ[end] = 1
ψ0 = vec(ψ)

Q = vec(z.^2 .- 2*z .+ 1)
Qmat = repeat(Q,1,N)
Zmat = repeat(2*(z .-1),1,N)
Ymat = repeat(1 ./y,1,N);

# Second Derivative
residuals = -0.5*( (Q.*(D2z*ψ) + 2*(z .-1).*(Dz*ψ) ).*Q/(4*R^2)
        + (Q/(2*R) ./y).*(Dz*ψ) )+ 0.5*K^2 *ψ ./y.^2 + ψ.^3 .- ψ
residuals[1] = 0; residuals[end] =0

while sum(abs.(residuals).^2) > 1e-12
    # Second Derivative
    residuals = -0.5*( (Q.*(D2z*ψ) + 2*(z .-1).*(Dz*ψ) ).*Q/(4*R^2)
            + (Q/(2*R) ./y).*(Dz*ψ) )+ 0.5*K^2 *ψ ./y.^2 + ψ.^3 .- ψ
    residuals[1] = 0; residuals[end] =0

    #println(sum(abs.(residuals).^2))


    Jacobi = -0.5*( (Qmat.*(D2z) + Zmat.*(Dz) ).*(Qmat/(4*R^2))
    + (Qmat/(2*R).*Ymat).*(Dz) ) + diagm(0 => 0.5*K^2 ./y.^2)+ diagm(0 => 3*ψ.^2 .- 1)


    Jacobi[1,:] = [1 zeros(1,N-1)]
    Jacobi[N,:] = [zeros(1,N-1) 1]


    Δ = vec(-4/7*(Jacobi\residuals))

    ψ = ψ + Δ
    ψ[1] = 0
    ψ[end] =1
end
    res = sum(abs.(residuals).^2)
    return y,ψ,res
end
