using VortexDistributions

# simple test field
Nx = 400; Ny = 400
Lx = 200; Ly = Lx
x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, Ny+1)[1:end-1];
psi0 = one.(x*y') |> complex

psi = Torus(copy(psi0), x, y);


xp = 0.0
yp = 0.0
vp = PointVortex(xp,yp,1)
sp = ScalarVortex(vp)

vortex!(psi, sp)

vfound = findvortices(psi)

vortm = keep_vortices(vfound)

@test length(vfound) == length(vortm)