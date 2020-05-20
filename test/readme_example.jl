## readme script
using Pkg;Pkg.activate("./")
using VortexDistributions, Plots
gr(xlabel="x",ylabel="y",transpose=true,legend=false)

## make a simple 2D test field
Nx = 400; Ny = Nx
Lx = 200; Ly = Lx
x = LinRange(-Lx / 2, Ly / 2, Nx); y = x
psi0 = one.(x*y') |> complex

## doubly periodic boundary conditions
psi = Torus(psi0,x,y)

# make a point vortex
pv = PointVortex(30.0,70.3,-1)

# make a scalar GPE vortex with exact core
spv = ScalarVortex(pv)
vortex!(psi,spv)

## make some more random vortices
vort = rand_vortex(10,psi)
vortex!(psi,vort)

using BenchmarkTools
@btime find_vortices(psi)
