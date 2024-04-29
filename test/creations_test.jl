using VortexDistributions

## Basic unit testing for coverage 
@test scalar_ansatz(1) == scalar_ansatz(-1)
f = Ansatz()
x = LinRange(0, 10, 100)
y = f.(x)
@test length(y) == length(x)
@test f.(3,4) == f(5)

y, ψ, res = VortexDistributions.gpecore_exact(1,2,100)
@test length(y) == 100

# create exact core
N = 100
L = 100
dx = L/N
x = LinRange(-L / 2, L / 2-dx, N); y = x
psi0 = one.(x*y') |> complex
psi = Torus(psi0,x,y)

# make a point vortex
pv = PointVortex(10.0,10.0,1)
nv = PointVortex(-10.0,-10.0,-1)
# make a scalar GPE vortex with exact core
spv = ScalarVortex(pv)
snv = ScalarVortex(nv)
vortex!(psi,spv)
vortex!(psi,snv)