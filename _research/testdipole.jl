using VortexDistributions
using Test, Plots

import VortexDistributions.thetad

## Periodic dipole phase
# Billam et al, PRL 112, 145301 (2014), Supplemental
# Summary:
# H using >= seems to agree with Billam. presumably atan handling of branch.
# paragraph 1 page 2 has a typo: shifted phase should involve x/(2π), y/(2π),
# as implemented in last line of thetad
H(x) = x >= 0. ? 1.0 : 0.0  # change to >=
shift(x,y) = x - y
tans(x,xk) = tan((shift(x,xk) - π)*0.5)
tanhs(x,xk,j) = tanh((shift(x,xk) + 2*π*j)*0.5)

function kernel(x,y,xp,yp,xn,yn,n)
    return atan(tanhs(y,yn,n)*tans(x,xn)) -
    atan(tanhs(y,yp,n)*tans(x,xp))
end

# domain [0,2*pi]; Eq (13); extra/(2π) in last line.
function thetad(x,y,xp,yp,xn,yn)
    s = 0.0
    for j ∈ -5:5
        s += kernel(x,y,xp,yp,xn,yn,j)
    end
    s += π*(H(shift(x,xp)) - H(shift(x,xn))) - y*(xp - xn)/(2*π)
    return s - x*H(abs(yp - yn) - π)/(2*pi) + y*H(abs(xp - xn) - π)/(2*pi) 
end




## distant dipole 1
    vp = PointVortex(-60,-60,1)
    vn = PointVortex(-60,60,-1)
    psi2 = Torus(copy(psi0),x,y)
    dip2 = ScalarVortex([vp,vn])
    periodic_dipole!(psi2,dip2)
    heatmap(x,y,abs2.(transpose(psi2.ψ)))
    heatmap(x,y,angle.(transpose(psi2.ψ)))


## distant dipole 2
    vp = PointVortex(-60.1,-60,1)
    vn = PointVortex(60,-60,-1)
    psi2 = Torus(copy(psi0),x,y)
    dip2 = ScalarVortex([vp,vn])
    periodic_dipole!(psi2,dip2)
    heatmap(x,y,abs2.(transpose(psi2.ψ)))
    heatmap(x,y,angle.(transpose(psi2.ψ)))


## distant dipole 3
    vp = PointVortex(-61,-60,1)
    vn = PointVortex(60,60,-1)
    psi2 = Torus(copy(psi0),x,y)
    dip2 = ScalarVortex([vp,vn])
    periodic_dipole!(psi2,dip2)
    heatmap(x,y,abs2.(transpose(psi2.ψ)))
    heatmap(x,y,angle.(transpose(psi2.ψ)))


## close dipole (works!)
    vp = PointVortex(20,10,1)
    vn = PointVortex(-10,.7,-1)
    dip = Dipole(vp,vn)
    @test typeof(dip) <: VortexGroup
    # simple test field
    Nx = 400; Ny = 400
    Lx = 200; Ly = Lx
    x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; 
    y = LinRange(-Ly/2,Ly/2, Ny+1)[1:end-1];
    psi0 = one.(x*y') |> complex
    psi = Torus(copy(psi0),x,y)
    dip2 = ScalarVortex([vp,vn])
    periodic_dipole!(psi,dip2)
    heatmap(x,y,abs2.(psi.ψ))
    heatmap(x,y,angle.(transpose(psi.ψ)))

## close dipole
    vp = PointVortex(-6.1,-6,1)
    vn = PointVortex(6,6,-1)
    psi2 = Torus(copy(psi0),x,y)
    dip2 = ScalarVortex([vp,vn])
    periodic_dipole!(psi2,dip2)
    heatmap(x,y,abs2.(psi2.ψ))
    heatmap(x,y,angle.(transpose(psi2.ψ)))


## large dipole. branch cut?
    vp = PointVortex(80,70,1)
    vn = PointVortex(-40,-30,-1)
    dip = Dipole(vp,vn)
    @test typeof(dip) <: VortexGroup
    # simple test field
    Nx = 400; Ny = 400
    Lx = 200; Ly = Lx
    x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]
    y = LinRange(-Ly/2,Ly/2, Ny+1)[1:end-1]
    psi0 = one.(x*y') |> complex
    psi = Torus(copy(psi0),x,y)
    dip2 = ScalarVortex([vp,vn])
    periodic_dipole!(psi,dip2)
    heatmap(x,y,abs2.(transpose(psi.ψ)))
    heatmap(x,y,angle.(transpose(psi.ψ)))


## summary:
# typo in Eq