

# unitary?
v0 = [.2 .4 1]
w0 = PointVortex(v0)
@test RawData(w0) == v0

v1 = [.2 .4 1;0.7 1.5 -1;-.3 1.2 1]
w1 = PointVortex(v1)
@test RawData(w1) == v1

w3 = randPointVortex(1000)
v3 = RawData(w3)
@test v3 == RawData(w3)

# single vortex creation and detection
@test foundNear(1)
@test foundNear(10)
@test foundNear(500)
