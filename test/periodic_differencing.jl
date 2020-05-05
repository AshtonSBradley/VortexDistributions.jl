# Lazy phase differencing
using Pkg;Pkg.activate(".")
using Test, Plots, BenchmarkTools, ShiftedArrays, Revise, VortexDistributions
gr(xlabel="x",ylabel="y",transpose=true,legend=false)

## unwrap make test phase
# test data
n=20
Nx = 400; Ny = Nx
Lx = 200; Ly = Lx
x = LinRange(-Lx / 2, Ly / 2, Nx); y = x
psi0 = one.(x*y') |> complex

# doubly periodic boundary conditions
psi = Torus(psi0,x,y)

# make a point vortex
pv = PointVortex(30.0,70.3,-1)

# make a scalar GPE vortex with exact core
spv = ScalarVortex(pv)
vortex!(psi,spv)

# make some more random vortices
vort = randVortex(n,psi)
vortex!(psi,vort)

ph = angle.(psi.ψ)

s1,s2 = size(ph)
ilag = mod1.(collect(1:s1) .-1, s1)

## unwrap test equivalence
testunwrap = unwrap(ph,1)
unwrap!(ph2,ph,1)
@test ph2 == testunwrap

## benchmark dim=1
@btime phu = unwrap(ph,1)
ph2 = deepcopy(ph)
@btime unwrap!(ph2,ph,1)

## benchmark dim=2
@btime phu = unwrap(ph,2)
ph2 = deepcopy(ph)
@btime unwrap!(ph2,ph,2)
