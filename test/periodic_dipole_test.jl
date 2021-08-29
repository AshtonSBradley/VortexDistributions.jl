using VortexDistributions

# simple test field
Nx = 400; Ny = 400
Lx = 200; Ly = Lx
x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, Ny+1)[1:end-1];
psi0 = one.(x*y') |> complex
psi1 = Torus(copy(psi0), x, y)
psi2 = Torus(copy(psi0), x, y)


# WEIRD BEHAVIOUR WHEN VORTICES ARE THIS CLOSE TO EACH OTHER, NEED TO INVESTIGATE 
# vp = PointVortex(.1,.3,1)
# vn = PointVortex(-.2,.7,-1)

xp = 1.1; yp = .3;
xn = -2.23; yn = .7;
vp = PointVortex(xp, yp, 1)
vn = PointVortex(xn, yn, -1)
dip = Dipole(vp,vn)
@test typeof(dip) <: VortexGroup

sp = ScalarVortex([vp, vn])
vortex!(psi1, sp)
vfound1 = findvortices(psi1)

periodic_dipole!(psi2, sp)
vfound2 = findvortices(psi2, periodic=true)

vpf1 = vfound1[2];
vnf1 = vfound1[1];
vpf2 = vfound2[2];
vnf2 = vfound2[1];

dx = Δ(psi1.x)
@test isapprox(vpf1.xv, xp, atol=dx/4)
@test isapprox(vpf1.yv, yp, atol=dx/4)
@test vpf1.qv == 1

@test isapprox(vnf1.xv, xn, atol=dx/4)
@test isapprox(vnf1.yv, yn, atol=dx/4)
@test vnf1.qv == -1

@test isapprox(vpf2.xv, xp, atol=dx/4)
@test isapprox(vpf2.yv, yp, atol=dx/4)
@test vpf2.qv == 1

@test isapprox(vnf2.xv, xn, atol=dx/4)
@test isapprox(vnf2.yv, yn, atol=dx/4)
@test vnf2.qv == -1


