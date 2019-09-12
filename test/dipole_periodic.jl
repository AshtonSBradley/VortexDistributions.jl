# doubly periodic dipole
# https://dx.doi.org/10.1103/PhysRevLett.112.145301
# should be eigenstate of Hamiltonian evolution

using Test, Plots, Revise, VortexDistributions


# Heaviside step
H(x) = x > 0. ? 1.0 : 0.0

# Translate to xi
T(x,xi) = x - xi

tanshift(x,xk) = tan((T(x,xk) - π)*0.5)
tanhshift(x,xk,j) = tanh((T(x,xk) + 2*π*j)*0.5)

# Vortex image kernel for image j
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
    for j = -5:5
        s += K(x,y,xp,yp,xn,yn,j)
    end
    return s - y*(xp - xn)/(2*π)
end

# test dipole
xp = 0.2
yp = .5
vp = PointVortex(xp,yp,1)

xn = 0.2
yn = -.5
vn = PointVortex(xn,yn,-1)

dip = [vp;vn]
θd(.1,.2,dip)

testphase = θd.(x,y',[dip])
testphase[1,:]

# test double periodicity
@test testphase[end,:] ≈ testphase[1,:]
@test testphase[:,1] ≈ testphase[:,end]

# test phase on "unit" length domain
x = LinRange(-pi,pi,300)
y = x
heatmap(x,y,θd.(x,y',[dip]),transpose=true)



# make ansatz wavefunction
# Default healing length ξ = 1.0
psi = one.(x*y') |> complex
dips = ScalarVortex(dip)
@. psi = psi * exp(im*θd(x,y',[dip]))
heatmap(x,y,angle.(psi),transpose=true)
