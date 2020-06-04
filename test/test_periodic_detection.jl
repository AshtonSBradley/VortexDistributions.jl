using Pkg;Pkg.activate(".")
using Test, Plots, Revise, VortexDistributions
gr(transpose=true)
# @test found_near(10)

## field with 2 points per healing length, test asymmetric
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

## detect
vort = find_vortices(psi)

## benchmark
using BenchmarkTools

## timing
@btime vort = find_vortices(psi)


## make a vortex at edge of grid
psi = Torus(copy(psi0),x,y)

xp = 100
yp = -10
vp = PointVortex(xp,yp,1)

xn = 100
yn = 10
vn = PointVortex(xn,yn,-1)

dip = ScalarVortex([vp;vn])

periodic_dipole!(psi,dip)


heatmap(angle.(psi.ψ))

## detect
vort = find_vortices(psi)
