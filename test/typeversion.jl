# using Parameters, Interpolations, Plots, Test
# using LinearAlgebra, ToeplitzMatrices, FileIO, JLD2
#
# abstract type Field end
# abstract type Vortex end
# abstract type CoreShape end
#
# @with_kw mutable struct Torus <: Field
#     ψ::Array{Complex{Float64},2}
#     x::Vector{Float64}
#     y::Vector{Float64}
# end
#
# @with_kw mutable struct Sphere <: Field
#     ψ::Array{Complex{Float64},2}
#     x::Vector{Float64}
#     y::Vector{Float64}
# end
#
# @with_kw mutable struct PointVortex <: Vortex
#     xv::Float64
#     yv::Float64
#     qv::Int64
# end
#
# PointVortex(v::Array{Float64,2}) = PointVortex.(v[:,1],v[:,2],v[:,3])
# RawData(v::PointVortex) = [v.xv v.yv v.qv]
# RawData(v::Array{PointVortex,1}) = reduce(vcat,RawData.(v))
#
# function uniform(a,b)
#     @assert a<b
#     return a + (b-a)*rand()
# end
#
# function uniform(a,b,n)
#     return uniform.(a*ones(n),b)
# end
#
# uniform() = uniform(-.5,.5)
# uniform(n) = uniform(-.5,.5,n)
#
# randcharge() = rand([-1 1],)
# randcharge(n) = rand([-1 1],n)
# randPointVortex(n) = PointVortex.(uniform(n), uniform(n), randcharge(n))
# randPointVortex() = PointVortex(uniform(),uniform(),randcharge())
#
# function randPointVortex(n,psi::Field)
#     @unpack ψ,x,y = psi
#     dx = x[2]-x[1]; dy = y[2]-y[1]
#     xi,xf = x[1]+2dx,x[end]-2dx
#     yi,yf = y[1]+2dy,y[end]-2dy
#     return PointVortex.(uniform(xi,xf,n),uniform(yi,yf,n),randcharge(n))
# end
#
# randPointVortex(psi::Field) = randPointVortex((),psi::Field)
#
# # Detection
# include("../src/detection.jl")
#
# function findwhere(A)
#     I = findall(!iszero,A)
#     v = A[I]
#     ix = [I[i][1] for i in eachindex(I)]
#     iy = [I[i][2] for i in eachindex(I)]
#     return ix,iy,v
# end
#
# function findvortices_jumps(psi::Field;shift=true)
#     @unpack x,y,ψ = psi
#     phase = angle.(ψ)
#
#     # Compute all nearest neighbour jumps
#     diffx = phasejumps(phase,1); diffy = phasejumps(phase,2)
#
#     # Compute all plaquette loop integrals with in-place memory recycling
#     circshift!(phase,diffx,(0,1))
#     diffx .-= phase; diffx .-= diffy
#     circshift!(phase,diffy,(1,0))
#     diffx .+= phase
#
#     # Find all windings
#     ixp,iyp,vp = findwhere(diffx .> 0.0)
#     xp = x[ixp]; yp = y[iyp]; np = length(vp)
#     ixn,iyn,vn = findwhere(diffx .< 0.0)
#     xn = x[ixn]; yn = y[iyn]; nn = length(vn)
#
#     if shift
#         dx = x[2]-x[1]; dy = y[2] - y[1]
#         xp .-= dx/2; yp .-= dy/2; xn .-= dx/2; yn .-= dy/2
#     end
#
#     vortices = [xn yn -vn; xp yp vp]
#
#     return PointVortex(vortices)
# end
#
# function findvortices_grid(psi::Torus;shift=true)
#     vort = findvortices_jumps(psi,shift=shift)
#     rawvort = RawData(vort)
#     vort = sortslices(rawvort,dims=1)
#     return PointVortex(vort)
# end
#
# function findvortices_grid(psi::Sphere;shift=true)
#     @unpack ψ = psi
#     windvals = phasejumps(angle.(ψ),2)
#     vort = findvortices_jumps(psi,shift=shift)
#     rawvort = RawData(vort)
#
#     # North pole winding. - sign means polar vortex co-rotating with w > 0
#     w1 = sum(windvals[1,:])
#     (sign(w1) != 0) && (rawvort = [rawvort; 0.0 0.0 -w1])
#     # South pole winding. + sign means co-rotating with w > 0
#     w2 = sum(windvals[end,:])
#     (sign(w2) != 0) && (rawvort = [rawvort; pi 0.0 w2])
#
#     vort = sortslices(vort,dims=1)
#     return PointVortex(vort)
# end
#
# function findvortices_interp(psi::Field)
#     vort = findvortices_grid(psi,shift=true)
#     vort = removeedgevortices(vort,psi)
#     #TODO: allow for interp with periodic data (here edges are stripped)
#
#     for (j,vortex) in enumerate(vort)
#         try
#         vortz,psiz = corezoom(vortex,psi)
#         vortz,psiz = corezoom(vortz,psiz)
#         vortz,psiz = corezoom(vortz,psiz)
#         vort[j] = vortz
#         catch nothing
#         end
#         return vort
#     end
# end
#
# """
#     vortices = findvortices(psi::T,interp=true) where T<:Field
#
# Locates vortices as 2π phase windings around plaquettes on a cartesian spatial grid. Uses an optimized plaquette method followed by recursive interpolation.
#
# Requires a 2D wavefunction ψ(x,y) on a cartesian grid specified by vectors x, y.
#
# `vortices` - array of vortex coordinates `xv,yv` and charges `qv`. Each row is of the form `[xv, yv, cv]`, and the array is sorted into lexical order according to the `xv` coordinates
# """
# function findvortices(psi::Field,interp=true)
#     if interp
#         return findvortices_interp(psi)
#     else
#         return findvortices_grid(psi)
#     end
# end
#
#
# """
# `vortices = removeedgevortices(vort::Array{PointVortex,1},x,y,edge=1)`
#
# Strips edgevortices due to periodic phase differencing."""
# function removeedgevortices(vort::Array{PointVortex,1},psi::Field,edge=1)
#     @unpack x,y = psi; dx=x[2]-x[1]; dy=y[2]-y[1]
#     keep = []
#     for j = 1:length(vort)
#         xi,yi,qi = RawData(vort[j])
#         xedge = isapprox(xi,x[1],atol=edge*dx) || isapprox(xi,x[end],atol=edge*dx)
#         yedge = isapprox(yi,y[1],atol=edge*dy) || isapprox(yi,y[end],atol=edge*dy)
#         not_edge = !(xedge || yedge)
#         not_edge && push!(keep,j)
#     end
#     return vort[keep]
# end
#
# """
#     vortz,psiz = corezoom(vortex,psi::T,winhalf=2,Nz=30) where T <: Field
#
# Uses local interpolation to resolve core location.
# """
# function corezoom(vortex::PointVortex,psi::T,winhalf=2,Nz=30) where T<:Field
#     @unpack ψ,x,y = psi
#     xv,yv,qv = RawData(vortex)
#     dx=x[2]-x[1]; dy=y[2]-y[1]
#     ixv = isapprox.(x,xv,atol=dx) |> findlast
#     iyv = isapprox.(y,yv,atol=dy) |> findlast
#     ixwin = (ixv-winhalf):(ixv+winhalf-1)
#     iywin = (iyv-winhalf):(iyv+winhalf-1)
#     xw = x[ixwin]; yw = y[iywin]; psiw = ψ[ixwin,iywin]
#     xz = LinRange(xw[1],xw[end],Nz)
#     yz = LinRange(yw[1],yw[end],Nz)
#     knots = (xw,yw)
#     itp = interpolate(knots, psiw, Gridded(Linear()))
#     psiz = itp(xz,yz)
#     ψv = T(psiz,xz |> Vector,yz |> Vector)
#     vortz = findvortices_grid(ψv,shift=false)
#     vortz = removeedgevortices(vortz,ψv)[1]
#     return vortz,ψv
# end
#
# corezoom(vortex::Array{PointVortex,1},psi::Field,winhalf=2,Nz=30) = corezoom(vortex[1],psi,2,30)

# Masking
# function findvortexmask(psi::T,R) where T<: Field
#     @unpack x, y, ψ = psi
#     ψ = circmask(psi, 1.1 * R)
#     ϕ = T(x, y, ψ)
#     vortices = findvortices(ϕ)
#
#     # TODO more tests
#     # remove vortices found outside mask boundary
#         for (i, xv) in enumerate(vortices[:,1]), yv in vortices[i,2]
#             (norm([xv,yv]) > R) && (vortices[i,:] = [0. 0. 0.])
#         end
#
#         vortfound = Array{Float64,2}[]
#         for (i, sv) in enumerate(vortices[:,3])
#             (sv != 0) && push!(vortfound, vortices[i,:]')
#         end
#     # currently will fail for more that one vortex (?)
#     return vortfound[1]
# end

# Helpers
# linspace(a,b,n) = LinRange(a,b,n) |> collect
#
# function edgemask!(psi::T) where T<:Field
#     @unpack ψ,x,y = psi
#     ψ[:,1] = zero(x) |> complex
#     ψ[:,end] = zero(x) |> complex
#     ψ[1,:] = zero(y') |> complex
#     ψ[end,:] = zero(y') |> complex
#     return T(ψ,x,y)
# end
#
# function circmask(psi::T,R) where T<:Field
#     @unpack ψ,x,y = psi
#     for j in eachindex(x), k in eachindex(y)
#         (x[j]^2+y[k]^2 > R^2) && (ψ[j,k] = complex(0.))
#     end
#     return T(x,y,ψ)
# end

# function circmask!(phi::T,psi::T,R) where T<: Field
#     @unpack ψ,x,y = psi
#     phi.ψ .= psi.ψ
#     for j in eachindex(x), k in eachindex(y)
#             (x[j]^2+y[k]^2 > R^2) && (phi.ψ[j,k] = complex(0.))
#     end
#     return nothing
# end
#
# function isinterior(a,b,x,y)
#     Lx = x[end]-x[1]; Ly = y[end] - y[1]
#     dx = x[2] - x[1]; dy = y[2] - y[1]
#     return (-Lx/2 + 2dx < a < Lx/2 - 2dx && -Ly/2 + 2dy < b < Ly/2 - 2dy)
# end

# #TODO update this function
# function checkvortexlocations(testvort,vortices,x,y,Nv)
# #check detection to 2 x grid resolution
#     dx = x[2] - x[1]; dy = y[2] - y[1]
#     vortfound = 0
#     for j=1:Nv
#         if (isapprox(testvort[j,1],vortices[j,1],atol=2dx) && isapprox(testvort[j,2],vortices[j,2],atol=2dy))
#         vortfound+=1
#         end
#     end
#     return vortfound
# end

# Creation
# include("../src/creation.jl")
# scalaransatz(x) = sqrt(x^2/(1+x^2))

using Parameters, Interpolations, Plots, Test
using LinearAlgebra, ToeplitzMatrices, FileIO, JLD2

include("../src/types.jl")
include("../src/pointvortex.jl")
include("../src/detection.jl")
include("../src/creation.jl")

# y,ψ,res = gpecore_exact(1)
# ψa = interpolate(tuple(y[1:99]), scalaransatz.(y[1:99]), Gridded(Linear()))
#
# ψi = interpolate(tuple(y[1:99]), ψ[1:99], Gridded(Linear()))
# @save "./src/ansatzcore.jld2" ψa
# @save "./src/exactcore.jld2" ψi

# @load "./src/exactcore.jld2" ψi
# @load "./src/ansatzcore.jld2" ψa


# r(x,y) = sqrt(x^2+y^2)
#
# struct Ansatz <: CoreShape
#     f::Interpolations.GriddedInterpolation{Float64,1,Float64,Gridded{Linear},Tuple{Array{Float64,1}}}
#     ξ::Float64
#     Λ::Float64
# end
#
# Ansatz() = Ansatz(ψa,1.0,0.8249)
#
# function (core::Ansatz)(x)
#     @unpack f,ξ,Λ = core
#     return f(Λ*x/ξ)
# end
# (core::Ansatz)(x,y) = core(r(x,y))
#
# struct Exact <: CoreShape
#     f::Interpolations.GriddedInterpolation{Float64,1,Float64,Gridded{Linear},Tuple{Array{Float64,1}}}
#     ξ::Float64
# end
# Exact() = Exact(ψi,1.0)
#
# function (core::Exact)(x)
#     @unpack ξ,f = core
#     return f(x/ξ)
# end
# (core::Exact)(x,y) = core(r(x,y))
#
# struct ScalarVortex{T<:CoreShape} <: Vortex
#     core::T
#     vort::PointVortex
# end
# ScalarVortex(vort::PointVortex) = ScalarVortex(Exact(),vort)
#
# function (s::ScalarVortex{T})(x,y) where T <: CoreShape
#     @unpack xv,yv,qv = s.vort
#     return s.core(x - xv, y - yv)*exp(im*qv*atan(y - yv,x - xv))
# end
#
# randPointVortex() = PointVortex(uniform(),uniform(),randcharge())
# randPointVortex(n) = PointVortex.(uniform(n), uniform(n), randcharge(n))
# function randPointVortex(psi::Field)
#     @unpack x,y = psi; a,b = first(x),last(x); c,d = first(y), last(y)
#     return PointVortex(uniform(a,b,),uniform(c,d,),randcharge())
# end
# function randPointVortex(n,psi::Field)
#     @unpack x,y = psi;
#     dx = x[2]-x[1]; dy = y[2] - y[1]
#     a,b = first(x)+2dx,last(x)-2dx; c,d = first(y)+2dy, last(y)-2dy
#     return PointVortex.(uniform(a,b,n),uniform(c,d,n),randcharge(n))
# end
#
# randScalarVortex() = ScalarVortex(randPointVortex())
# randScalarVortex(n) = ScalarVortex.(randPointVortex(n))
# randScalarVortex(psi::Field) = ScalarVortex(Exact(),randPointVortex(psi))
# randScalarVortex(n,psi::Field) = ScalarVortex.([Exact()],randPointVortex(n,psi))
#
# randVortex() = randScalarVortex()
# randVortex(n) = randScalarVortex(n)
# randVortex(n,psi::Field) = randScalarVortex(n,psi)
#
# function vortex!(psi::F,vort::ScalarVortex{T}) where {T <: CoreShape, F<:Field}
#     @unpack ψ,x,y = psi
#     @. ψ *= vort.(x,y')
#     @pack! psi = ψ
# end
#
# function vortex!(psi::F,vort::Array{S}) where {F <: Field, S <: Vortex}
#     for j in eachindex(vort)
#         vortex!(psi,vort[j])
#     end
# end

#=
# tests
x = LinRange(0,1,500)
ansatz = Ansatz(ψa,1.,0.8249)
ansatz(.1)
@time ansatz.(x)
@time ansatz.(x,x')
@time Ansatz(ψa,1.2,0.8249).(x,x')

exact = Exact()

exact(.1)
@time exact.(x)
@time exact.(x,x')
@time Exact(ψi,1.2).(x,x')
@time ψi.(x)
@time scalaransatz.(x)

vort1 = randPointVortex()
scalar_exact = ScalarVortex(exact,vort1)
scalar_ansatz = ScalarVortex(ansatz,vort1)
scalar_test = ScalarVortex(vort1)

scalar_ansatz.(x,x')
scalar_exact.(x,x')
scalar_test.(x,x')
@test ScalarVortex(vort1).(x,x') == scalar_test.(x,x')

using Plots
gr(transpose=true,aspectratio=1)
N = 200
L = 100.0
x = LinRange(-L/2,L/2,N)
y = x

psi0 = one.(x*y') |> complex; psi = Torus(psi0,x,y)
vtest = randPointVortex(psi)
vort1 = ScalarVortex(exact,vtest)
vort2 = ScalarVortex(Exact(),vtest)

vort3 = ScalarVortex(vtest)
psi0 = one.(x*y') |> complex; psi = Torus(psi0,x,y)
vortex!(psi,vort3)
@unpack ψ = psi
heatmap(x,y,angle.(ψ))
heatmap(x,y,abs2.(ψ))

vfound = findvortices(psi)

vort4 = ScalarVortex(Ansatz(),vtest)
psi0 = one.(x*y') |> complex
psi = Torus(psi0,x,y)
vortex!(psi,vort4)
@unpack ψ = psi
heatmap(x,y,angle.(ψ))
heatmap(x,y,abs2.(ψ))

# timing
x = LinRange(0, 6, 300)

vort4.(x,x')
@time vort4.(x,x')
@time vort3.(x,x')
@time vort2.(x,x')



# ==== old tests



# exact core
# y,ψ,res = gpecore_exact(1)
# @load "./src/exactcore.jld2" ψi

# Testing
gr(transpose=true,aspectratio=1)
N = 200
L = 100.0
x = LinRange(-L/2,L/2,N)
y = x

psi0 = one.(x*y') |> complex
psi = Torus(psi0,x,y)
vtest = randVortex(psi)

@time vortex!(psi,vtest)
@unpack ψ = psi
heatmap(x,y,angle.(ψ))
heatmap(x,y,abs2.(ψ))

# Testing
gr(transpose=true,aspectratio=1)
N = 200
L = 100.0
x = LinRange(-L/2,L/2,N)
y = x

psi0 = one.(x*y') |> complex
psi = Torus(psi0,x,y)
vtest = randVortex(psi)
vortex!(psi,vtest)
@unpack ψ = psi
heatmap(x,y,angle.(ψ))
heatmap(x,y,abs2.(ψ))

# test detection
vgridraw = findvortices_grid(psi,shift=false)
vgridraw = removeedgevortices(vgridraw,psi)
scatter!([RawData(vgridraw)[1]],[RawData(vgridraw)[2]])

vgrid = findvortices_grid(psi)
vgrid = removeedgevortices(vgrid,psi)
scatter!([RawData(vgrid)[1]],[RawData(vgrid)[2]])
=#
