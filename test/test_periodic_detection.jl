using Pkg;Pkg.activate(".")
using Test, Plots, BenchmarkTools, Revise, VortexDistributions
gr(xlabel="x",ylabel="y",transpose=true,legend=false)

## readme benchmark
# make a periodic dipole at boundary
Nx = 400; Ny = Nx
Lx = 200; Ly = Lx
x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = x
psi0 = one.(x*y') |> complex

# doubly periodic boundary conditions
psi = Torus(psi0,x,y)

# test dipole
xp = minimum(y)
yp = 0.1
vp = PointVortex(xp,yp,1)

xn = minimum(y)
yn = 3*pi
vn = PointVortex(xn,yn,-1)

dip = ScalarVortex([vp;vn])

periodic_dipole!(psi,dip)

heatmap(x,y,angle.(psi.ψ))

vort = findvortices(psi,interp=false)

vort = findvortices(psi)
# NOTE: why reverted to grid only? dx/2,dy/2 offset
# TODO: a systematic error to be tracked down

## detect
@btime vort = findvortices(psi,interp=false)

## run benchmark (cf 4ms)

@btime vort = findvortices(psi)

vort = findvortices(psi)
rawData(vort)











## test set
# unitary?
v0 = [.2 .4 1]
w0 = PointVortex(v0)
@test rawData(w0) == v0

v1 = [.2 .4 1;0.7 1.5 -1;-.3 1.2 1]
w1 = PointVortex(v1)
@test rawData(w1) == v1

w3 = randPointVortex(1000)
v3 = rawData(w3)
@test v3 == rawData(w3)

# single vortex creation and detection
@test foundNear(1)
@test foundNear(30)
