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

xn = 3.24
yn = 10.21
vn = PointVortex(xn,yn,-1)

dip = ScalarVortex([vp;vn])
periodic_dipole!(psi,dip)

heatmap(x,y,angle.(psi.ψ))

vfound = findvortices(psi)
@show vfound

## best test: find a known vortex
# Periodic boundary conditions
psi = Torus(copy(psi0),x,y)

xp = 1.133
yp = 2.787
vp = PointVortex(xp,yp,1)

sp = ScalarVortex(vp)

vortex!(psi,sp)
vfound = findvortices(psi);@show vortex_array(vfound)

heatmap(x,y,angle.(psi.ψ))
## test zoom_grid

psiw,xw,yw = zoom_grid(psi.ψ,x,y,vortex_array(vp)[1:2]...)
heatmap(xw,yw,angle.(psiw))


## readme example
@btime findvortices(psi)
