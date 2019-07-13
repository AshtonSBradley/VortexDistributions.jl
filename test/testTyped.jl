# test typed version
using Test

include("typeversion.jl")

# unitary?
v0 = [.2 .4 1]
w0 = PointVortex(v0)
@test RawVortex(w0) == v0

v1 = [.2 .4 1;0.7 1.5 -1;-.3 1.2 1]
w1 = PointVortex(v1)
@test RawVortex(w1) == v1


w3 = randvortex(1000)
v3 = RawVortex(w3)
@test PointVortex(v3) == w3

# fast. really, really fast!
@time v3 = RawVortex(w3)
@time w3 = PointVortex(v3)
