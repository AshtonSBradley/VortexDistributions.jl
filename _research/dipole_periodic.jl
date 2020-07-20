## doubly periodic dipole
# https://dx.doi.org/10.1103/PhysRevLett.112.145301
# should be eigenstate of Hamiltonian evolution

using Pkg;Pkg.activate(".")
using Test, Plots, VortexDistributions

## test dipole
xp = 20
yp = 3
vp = PointVortex(xp,yp,1)

xn = 1.1
yn = -3
vn = PointVortex(xn,yn,-1)



Thetad(.12,.11,xp,yp,xn,yn)

x = LinRange(-π,π,300)
y = x

#TODO add missing method for Thetad(x,y,dipole::Array{PointVortex})

## test double periodicity
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


#--- add new method for making periodic dipoles
using Test, Plots, Revise, VortexDistributions

x = LinRange(-50,50,300)
y = x

# test dipole
xp = 20
yp = 3
vp = PointVortex(xp,yp,1)

xn = 1.1
yn = -3
vn = PointVortex(xn,yn,-1)

dip = ScalarVortex([vp;vn])

function periodic_dipole!(psi::F,dip::Array{ScalarVortex{T},1}) where {T <: CoreShape, F<:Field}
    @assert length(dip) == 2
    @assert dip[1].vort.qv + dip[2].vort.qv == 0
    @assert hypot(dip[1].vort.xv-dip[2].vort.xv,dip[1].vort.yv-dip[2].vort.yv) >= 2*pi
    @unpack ψ,x,y = psi
    (dip[1].vort.qv > 0) ? (jp = 1;jn = 2) : (jp = 2;jn = 1)
    vp = rawData(dip[jp].vort)[1:2]
    vn = rawData(dip[jn].vort)[1:2]
    @. ψ *= abs(dip[jn](x,y')*dip[jp](x,y'))
    ψ .*= exp.(im*Thetad(x,y,vp...,vn...))
    @pack! psi = ψ
end

psivort = Torus(one.(x.*y'),x,y)

heatmap(abs2.(psivort.ψ))

periodic_dipole!(psivort,dip)

heatmap(x,y,abs2.(psivort.ψ))

heatmap(x,y,angle.(psivort.ψ))

#check that phase is periodic (note we don't enforce periodic density)

testpsi = angle.(psivort.ψ)
@test testpsi[:,1] ≈ testpsi[:,end]
@test testpsi[1,:] ≈ testpsi[end,:]
