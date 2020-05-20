## activate
using Pkg;
Pkg.activate("./")
using Test, Revise, VortexDistributions

## unitary?
v0 = [.2 .4 1]
w0 = PointVortex(v0)
@test vortex_array(w0) == v0

v1 = [.2 .4 1;0.7 1.5 -1;-.3 1.2 1]
w1 = PointVortex(v1)
@test vortex_array(w1) == v1

w3 = rand_pointvortex(1000)
v3 = vortex_array(w3)
@test v3 == vortex_array(w3)

## single vortex creation and detection
@test found_near(1)
@test found_near(30)

## create a dipole
vp = PointVortex(.1,.3,1)
vn = PointVortex(-.2,.7,-1)
dip = Dipole(vp,vn)

## Test positions
@test xpos(vp) == 0.1
@test ypos(vp) == 0.3

vv = [vp; vn]
@test xpos(vv) == [.1;-.2]
@test ypos(vv) == [.3;.7]
pos(vp)
pos(vv)
charge(vv)

## create a cluster
# cl = Cluster()
