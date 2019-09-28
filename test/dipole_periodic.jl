# doubly periodic dipole
# https://dx.doi.org/10.1103/PhysRevLett.112.145301
# should be eigenstate of Hamiltonian evolution

using Test, Plots, Revise, VortexDistributions

# test dipole
xp = 0.2
yp = .2
vp = PointVortex(xp,yp,1)

xn = 1.1
yn = -.2
vn = PointVortex(xn,yn,-1)

dip = [vp;vn]

Thetad(.12,.11,xp,yp,xn,yn)

x = LinRange(-π,π,300)
y = x

#TODO add missing method for Thetad(x,y,dipole::Array{PointVortex})
# test double periodicity
testpsi = exp.(im*Thetad(x,y,xp,yp,xn,yn))
@test testpsi[end,:] ≈ testpsi[1,:]
@test testpsi[:,1] ≈ testpsi[:,end]
heatmap(x,y,angle.(testpsi),transpose=true)

# test domain
x = LinRange(-3π,3π,300)
y = LinRange(-4π,4π,300)

# test dipole
xp = 2.1*pi
yp = 2*pi

xn = -2*pi
yn = -2.4*pi

testphase3 = Thetad(x,y,xp,yp,xn,yn)
heatmap(x,y,testphase3,transpose=true)

psivort = Torus(exp.(im*testphase3),x,y)
vortices = findvortices(psivort)

vraw = rawData(vortices)

testpsi = exp.(im*testphase3)
@test testpsi[:,1] ≈ testpsi[:,end]
@test testpsi[1,:] ≈ testpsi[end,:]
