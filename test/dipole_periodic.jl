# doubly periodic dipole
# https://dx.doi.org/10.1103/PhysRevLett.112.145301
# should be eigenstate of Hamiltonian evolution

using Revise, VortexDistributions
# test dipole propagates along positive x
# vd = 1.

xp = 0.
xn = 0.
yp = 0.5
yn = -0.5

vp = PointVortex(xp,yp,1)
vn = PointVortex(xn,yn,-1)
dip = [vp;vn]

# default healing length ξ = 1.0
dips = ScalarVortex(dip)

# Heaviside step
H(x) = x > 0. ? 1.0 : 0.0

# Translate to xi
T(x,xi) = x - xi

tanshift(x,xk) = tan(T(x,xk) - π)
tanhshift(y,yk,j) = tanh((T(y,yk) + 2*π*j)*0.5)

# vortex image kernel for image j
function K(x,y,xp,yp,xn,yn,j)
    return atan(tanhshift(y,yn,j)*tanshift(x,xn)) -
    atan(tanhshift(y,yp,j)*tanshift(x,xp)) +
    π*(H(T(x,xp)) - H(T(x,xn)))
end

# dimensionless form
function θd(x,y,dip)
    vp,vn = dip
    xp,yp,_ = rawData(vp)
    xn,yn,_ = rawData(vn)
    s = 0.0
    for j = 1:5
        s += K(x,y,xp,yp,xn,yn,j)
    end
    return s - y*(xp - xn)/(2*π)
end

θd(.1,.1,dip)

# test wavefunction
using Plots
x = LinRange(-pi,pi,300)
y = x
psi = one(x*y') |> complex

θd.(x,y',[dip])

@. psi = psi * exp(im*θd(x,y',[dip]))

heatmap(x,y,angle.(psi))
