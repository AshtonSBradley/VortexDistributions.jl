# simple test field
using VortexDistributions, Test
Nx = 400; Ny = 400
Lx = 200; Ly = Lx
x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, Ny+1)[1:end-1];
psi = one.(x*y') |> complex

vfound = findvortices(Torus(psi,x,y), periodic=true)
