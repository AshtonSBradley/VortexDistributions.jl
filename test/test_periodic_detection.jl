using Pkg;Pkg.activate(".")
using Test, Plots, Revise, VortexDistributions

## readme test
# make a periodic dipole at boundary
Nx = 400; Ny = Nx
Lx = 200; Ly = Lx
x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = x
psi0 = one.(x*y') |> complex

dx = x[2]-x[1]
# doubly periodic boundary conditions
psi = Torus(psi0,x,y)

# test dipole
xp = 1.1
yp = 2.1
vp = PointVortex(xp,yp,1)

xn = 2.1
yn = 3.1
vn = PointVortex(xn,yn,-1)

dip = ScalarVortex([vp;vn])

periodic_dipole!(psi,dip)

heatmap(x,y,angle.(psi.ψ))

vort = findvortices(psi)
