using VortexDistributions
using Test

@testset "Point vortex " begin

    # unitary
    v0 = [.2 .4 1]
    w0 = PointVortex(v0)
    @test vortex_array(w0) == v0

    v1 = [.2 .4 1;0.7 1.5 -1;-.3 1.2 1]
    w1 = PointVortex(v1)
    @test vortex_array(w1) == v1

    w3 = rand_pointvortex(1000)
    v3 = vortex_array(w3)
    @test v3 == vortex_array(w3)

end

@testset "Single vortex accuracy" begin

    # simple test field
    Nx = 400; Ny = 400
    Lx = 200; Ly = Lx
    x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, Ny+1)[1:end-1];
    psi0 = one.(x*y') |> complex
    psi = Torus(copy(psi0),x,y)
    xp = 1.133
    yp = 2.787
    vp = PointVortex(xp,yp,1)
    sp = ScalarVortex(vp)
    vortex!(psi,sp)
    vfound = findvortices(psi)
    @test isapprox(vfound[1].xv,xp,atol=0.01)
    @test isapprox(vfound[1].yv,yp,atol=0.01)
end

@testset "Multivortex creation and detection" begin
    @test found_near(30)
end

@testset "Periodic dipole" begin
    vp = PointVortex(.1,.3,1)
    vn = PointVortex(-.2,.7,-1)
    dip = Dipole(vp,vn)
    @test typeof(dip) <: VortexGroup
end
