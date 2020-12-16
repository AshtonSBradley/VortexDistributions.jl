using VortexDistributions
using Test, Plots

## Billam et al, PRL 112, 145301 (2014), Supplemental
## close dipole
    vp = PointVortex(1,.3,1)
    vn = PointVortex(-1,.7,-1)
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
    heatmap(x,y,angle.(psi.ψ))

## distant dipole 1
    vp = PointVortex(-60,-60,1)
    vn = PointVortex(-60,60,-1)
    psi2 = Torus(copy(psi0),x,y)
    dip2 = ScalarVortex([vp,vn])
    periodic_dipole!(psi2,dip2)
    heatmap(x,y,abs2.(psi2.ψ))
    heatmap(x,y,angle.(psi2.ψ))

## distant dipole 2
    vp = PointVortex(-60,-60,1)
    vn = PointVortex(60,-60,-1)
    psi2 = Torus(copy(psi0),x,y)
    dip2 = ScalarVortex([vp,vn])
    periodic_dipole!(psi2,dip2)
    heatmap(x,y,abs2.(psi2.ψ))
    heatmap(x,y,angle.(psi2.ψ))