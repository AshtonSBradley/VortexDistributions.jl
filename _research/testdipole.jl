using VortexDistributions
using Test, Plots

    # simple test field
    Nx = 400; Ny = 400
    Lx = 200; Ly = Lx
    x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, Ny+1)[1:end-1];
    psi0 = one.(x*y') |> complex
    psi = Torus(copy(psi0),x,y)
    xp = 1.133
    yp = 2.787
    vp = PointVortex(xp,yp,1)
    sp = ScalarVortex(vp)

## distant dipole 1
    vp = PointVortex(-60,-60,1)
    vn = PointVortex(-60,60,-1)
    psi2 = Torus(copy(psi0),x,y)
    dip2 = ScalarVortex([vp,vn])
    periodic_dipole!(psi2,dip2)
    heatmap(x,y,abs2.(transpose(psi2.ψ)))
    heatmap(x,y,angle.(transpose(psi2.ψ)))

    vort = findvortices(psi2)

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