#find single vortex position to machine precision
using VortexDistributions, LinearAlgebra, Plots, Revise

#make a large homogeneous wavefunction (ξ=1)
Lx = 100.; Ly = 100.
Nx = 1000; Ny = 1000
x = linspace(-Lx/2,Lx/2,Nx+1); y = linspace(-Ly/2,Ly/2,Ny+1)
x = x[1:end-1];y = y[1:end-1]
dx = x[2]-x[1]; dy = y[2]-y[1]

# Vortex at origin
x1,y1 = (0.,0.)
vort = [x1 y1 1]

#construct vortex wavefunction (core ansatz)
psi = ones(size(x.*y')) |> complex

makeallvortices!(psi,vort,x,y,1.0)
heatmap(x,y,angle.(psi),xlabel="x",ylabel="y",transpose = true)

np,nn,nt,vortices = findvortices(psi,x,y)
vortices = remove_edgevortices(vortices,x,y)

# Random in the box interior
#place the vortices
x1,y1 = 0.8*Lx*(rand(2).-.5)
vort = [x1 y1 1]

#construct vortex wavefunction (core ansatz)
psi = ones(size(x.*y')) |> complex

makeallvortices!(psi,vort,x,y,1.0)
heatmap(x,y,angle.(psi),xlabel="x",ylabel="y",transpose = true)

np,nn,nt,vortices = findvortices(psi,x,y)
vortices = remove_edgevortices(vortices,x,y)

xv,yv = vortices[1:2]

x .== xv
ix = findfirst(x .== xv+dx/2)
x[ix]-dx/2
iy = findfirst(y .== yv+dy/2)
y[iy]
