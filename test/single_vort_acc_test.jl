using VortexDistributions, Test

# @testset "Single vortex accuracy" begin

# simple test field
Nx = 400; Ny = 400
Lx = 200; Ly = Lx
x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, Ny+1)[1:end-1];
psi0 = one.(x*y') |> complex

psi1 = Torus(copy(psi0),x,y)
psi2 = Torus(copy(psi0),x,y)
xp = 1.133
yp = 2.787
vp = PointVortex(xp,yp,1)
sp = ScalarVortex(vp)

vortex!(psi1,sp)
vfound1 = findvortices(psi1)
@test isapprox(vfound1[1].xv,xp,atol=0.01)
@test isapprox(vfound1[1].yv,yp,atol=0.01)


vortex!(psi2, vp)
vfound2 = findvortices(psi2)
@test isapprox(vfound2[1].xv,xp,atol=0.01)
@test isapprox(vfound2[1].yv,yp,atol=0.01)

@test (vfound1[1].xv === vfound2[1].xv && vfound1[1].yv === vfound2[1].yv && vfound1[1].qv === vfound2[1].qv)
# end
