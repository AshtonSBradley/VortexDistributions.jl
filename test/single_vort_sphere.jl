using VortexDistributions

# simple test field
Nx = 400; Ny = 400
Lx = 200; Ly = Lx
x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, Ny+1)[1:end-1];
psi0 = one.(x*y') |> complex

psi = Sphere(copy(psi0),x,y)

xp = 1.133
yp = 2.787
vp = PointVortex(xp,yp,1)
sp = ScalarVortex(vp)

vortex!(psi, sp)

vfound = findvortices(psi)

dx = Δ(psi.x)

@test isapprox(vfound[1].xv, xp, rtol=dx/4)
@test isapprox(vfound[1].yv, yp, rtol=dx/4)
@test vfound[1].qv == 1