using VortexDistributions, Test

# @testset "Periodic dipole" begin
vp = PointVortex(.1,.3,1)
vn = PointVortex(-.2,.7,-1)
dip = Dipole(vp,vn)
@test typeof(dip) <: VortexGroup
# end
