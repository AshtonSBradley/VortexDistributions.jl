# doubly periodic dipole
# https://dx.doi.org/10.1103/PhysRevLett.112.145301
# should be eigenstate of Hamiltonian evolution

using Test, Plots, Revise, VortexDistributions

# Methods
H(x) = x > 0. ? 1.0 : 0.0
T(x,xi) = x - xi
tans(x,xk) = tan((T(x,xk) - π)*0.5)
tanhs(x,xk,j) = tanh((T(x,xk) + 2*π*j)*0.5)

function K(x,y,xp,yp,xn,yn,j)
    return atan(tanhs(y,yn,j)*tans(x,xn)) -
    atan(tanhs(y,yp,j)*tans(x,xp))
end

# Dimensionless form
function θd(x,y,dip)
    vp,vn = dip
    xp,yp,_ = rawData(vp)
    xn,yn,_ = rawData(vn)
    s = 0.0
    for j = -5:5
        s += K(x,y,xp,yp,xn,yn,j)
    end
    return s + π*(H(T(x,xp)) - H(T(x,xn))) - y*(xp - xn)/(2*π)
end

# test dipole
xp = 1.5
yp = .2
vp = PointVortex(xp,yp,1)

xn = 1.7
yn = -.2
vn = PointVortex(xn,yn,-1)

dip = [vp;vn]
θd(.1,.2,dip)

x = LinRange(-π,π,300)
y = x

# test double periodicity
testphase = θd.(x,y',[dip])
@test testphase[end,:] ≈ testphase[1,:]
@test testphase[:,1] ≈ testphase[:,end]
heatmap(x,y,testphase,transpose=true)

# simplify the branch cut (only need phase up to + n*2*π)
testphase2 = angle.(exp.(im*θd.(x,y',[dip])))
@test testphase2[end,:] ≈ testphase2[1,:]
@test testphase2[:,1] ≈ testphase2[:,end]
heatmap(x,y,testphase2,transpose=true)



# Now for arbitrary domains and dipole sizes:
function thetad(x,y,xp,yp,xn,yn)
    s = 0.0
    for j = -5:5
        s += K(x,y,xp,yp,xn,yn,j)
    end
    s += π*(H(T(x,xp)) - H(T(x,xn))) - y*(xp - xn)/(2*π)
    return s - x*H(abs(yp - yn) - π) + y*H(abs(xp - xn) - π)
end

function Thetad(x,y,xp,yp,xn,yn)
    Lx = x[end]-x[1]
    Ly = y[end]-y[1]
    return @. angle(exp(im*thetad.(x*2*pi/Lx,y'*2*pi/Ly,xp*2*pi/Lx,yp*2*pi/Ly,xn*2*pi/Lx,yn*2*pi/Ly)))
end

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
