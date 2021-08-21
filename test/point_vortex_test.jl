using VortexDistributions, Test

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