using Pkg;Pkg.activate(".")
using Test, Plots, BenchmarkTools, Revise, VortexDistributions
gr(transpose=true)

## field with 2 points per healing length
Nx = 400; Ny = 400
Lx = 200; Ly = Lx
x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, Ny+1)[1:end-1];
psi0 = one.(x*y') |> complex

## make vortex dipole (cant use this for single vortex accuracy)
# Periodic boundary conditions
psi = Torus(copy(psi0),x,y)

xp = 1.7
yp = 2.1
vp = PointVortex(xp,yp,1)

xn = 3.25
yn = 10.
vn = PointVortex(xn,yn,-1)

dip = ScalarVortex([vp;vn])
periodic_dipole!(psi,dip)

heatmap(x,y,angle.(psi.ψ))

vfound = findvortices(psi)

## single vortex
# Periodic boundary conditions
psi = Torus(copy(psi0),x,y)

xp = 1.13
yp = 2.78
vp = PointVortex(xp,yp,1)

sp = ScalarVortex(vp)

vortex!(psi,sp)

heatmap(x,y,angle.(psi.ψ))
## test zoom_grid

psiw,xw,yw = zoom_grid(psi.ψ,x,y,vortex_array(vp)[1:2]...)
heatmap(xw,yw,angle.(psiw))

## test zoom_interp
# TODO mystery of why shift has no effect but to make x grid resolved...
## test findvortices

@btime vfound = findvortices(psi)
