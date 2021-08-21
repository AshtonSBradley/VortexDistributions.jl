using VortexDistributions, Test

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
