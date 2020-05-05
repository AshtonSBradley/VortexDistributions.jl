using Pkg;Pkg.activate(".")
using Test, Plots, Revise, VortexDistributions

@test foundNear(10)

## field with 2 points per healing length
Nx = 400; Ny = Nx
Lx = 200; Ly = Lx
x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = x
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

periodic_dipole!(psi,dip)

# heatmap(x,y,angle.(psi.ψ))

## detect
vort = findvortices(psi)

## benchmark
using BenchmarkTools

@btime vort = findvortices(psi)
