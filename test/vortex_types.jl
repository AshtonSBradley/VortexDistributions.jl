using VortexDistributions

# simple test field
Nx = 400; Ny = 400
Lx = 200; Ly = Lx
x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, Ny+1)[1:end-1];
psi0 = one.(x*y') |> complex
psi = Torus(copy(psi0),x,y)
xp = 1.133
yp = 2.787
vp = PointVortex(xp,yp,1)

## Test differnt calls to ScalarVortex()
sp = ScalarVortex(vp)
sp1 = ScalarVortex(1.0, vp)
@test sp1[1].core === sp.core
@test sp1[1].vort === sp.vort

sp2 = ScalarVortex(0.8, vp)
@test !(sp2[1].core===sp.core)
@test sp2[1].vort === sp.vort


sp3 = ScalarVortex([1.0, 1.0], [vp])
@test sp3[1].core === sp3[2].core
@test sp3[1].vort === sp3[2].vort

sp4 = ScalarVortex([1.0, 0.8], [vp])
@test !(sp4[1].core === sp4[2].core)
@test (sp4[1].vort === sp4[2].vort)

sp5 = ScalarVortex([vp])
@test sp5[1].core === sp1[1].core
@test sp5[1].vort === sp1[1].vort

rs1 = rand_scalarvortex()
@test typeof(rs1) <: ScalarVortex{Exact}

rs2 = rand_scalarvortex(10)
@test length(rs2) == 10
@test typeof(rs2[1]) <: ScalarVortex{Exact}

rs3 = rand_scalarvortex(10, psi)
@test typeof(rs3) == typeof(rs2)
@test length(rs3) == 10

rs4 = rand_scalarvortex(psi)
@test typeof(rs4) <: ScalarVortex{Exact}

rv1 = rand_vortex()
@test typeof(rv1) <: ScalarVortex{Exact}

rv2 = rand_vortex(10)
@test length(rv2) == 10
@test typeof(rv2[1]) <: ScalarVortex{Exact}

rv3 = rand_scalarvortex(10, psi)
@test typeof(rv3) == typeof(rv2)
@test length(rv3) == 10
