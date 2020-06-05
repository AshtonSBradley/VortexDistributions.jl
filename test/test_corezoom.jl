using Pkg;Pkg.activate(".")
using Test, Plots, Revise, VortexDistributions
gr(transpose=true)

## definition
using ParameterizedFunctions, Interpolations
"""
    vortz,psiz = corezoom(vortex,ψ<:Field,winhalf=2,Nz=30)

Uses local interpolation to resolve core location.
"""
function corezoom(vortex::PointVortex,psi::T,winhalf=2,Nz=30) where T<:Field
    @info "made it here..."
    @unpack ψ,x,y = psi
    xv,yv,qv = vortex.x,vortex.y.vortex.q
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
    vortz = find_vortices_grid(ψv,shift=true)
    vortz = remove_edge_vortices(vortz,ψv)[1]
    return vortz,ψv
end

## field with 2 points per healing length
Nx = 400; Ny = 400
Lx = 200; Ly = Lx
x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, Ny+1)[1:end-1];
psi0 = one.(x*y') |> complex

## make dipole
# Periodic boundary conditions
psi = Torus(copy(psi0),x,y)

xp = 1.7
yp = 2.1
vp = PointVortex(xp,yp,1)

xn = pi
yn = 10.
vn = PointVortex(xn,yn,-1)

dip = ScalarVortex([vp;vn])
# dip = Dipole(vp,vn)
# vortex_array(dip.vp)[1:2]

periodic_dipole!(psi,dip)

heatmap(x,y,angle.(psi.ψ))

## detect:  this must give sub-grid detection

vort = findvortices(psi)
@show vort[1]

## benchmark
using BenchmarkTools

## timing
@btime vort = findvortices(psi)
