using Pkg;Pkg.activate(".")
using Test, Plots, BenchmarkTools, Revise
using VortexDistributions
gr(transpose=true)

## field with 2 points per healing length
Nx = 400; Ny = 400
Lx = 200; Ly = Lx
x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, Ny+1)[1:end-1];
psi0 = one.(x*y') |> complex

## make dipole
# Periodic boundary conditions
psi = Torus(copy(psi0),x,y)

xp = 1.7
yp = 2.1
vp = PointVortex(xp,yp,1)

xn = pi
yn = 10.
vn = PointVortex(xn,yn,-1)

v1 = ScalarVortex(vp)

vortex!(psi,v1)

heatmap(x,y,angle.(psi.ψ))

## detect: grid limited?
vort = findvortices(psi)
@show vort[1]


## make a field with one boundary vortex
# Periodic boundary conditions
psi = Torus(copy(psi0),x,y)

xp = 100
yp = 0
vp = PointVortex(xp,yp,1)

xn = pi
yn = 10.
vn = PointVortex(xn,yn,-1)

v1 = ScalarVortex(vp)

vortex!(psi,v1)

heatmap(x,y,angle.(psi.ψ))

## periodic
import VortexDistributions:findvortices, zoom_interp, zoom_grid
function findvortices(psi::Field)
    @unpack ψ,x,y = psi
    vort = findvortices_grid(psi)
    # vort = remove_vortices_edge(vort,psi) #periodic: don't remove

    for (j,vortex) in enumerate(vort)
        v = try
        psi_int,xint,yint = zoom_interp(ψ,x,y,vortex.xv,vortex.yv) #periodic: peridic indices here
        v1 = findvortices_grid(Torus(psi_int,xint,yint))
        vint = remove_vortices_edge(v1,Torus(psi_int,xint,yint))[1]
        psi_int,xint,yint = zoom_interp(psi_int,xint,yint,vint.xv,vint.yv)
        v1 = findvortices_grid(Torus(psi_int,xint,yint))
        vint = remove_vortices_edge(v1,Torus(psi_int,xint,yint))[1]
        catch nothing
        end
        vort[j] = v     # TODO a fallback for not found?
        # periodic: map coordinates to fundamental box
    end
    return vort
end

"""
    psiz,xz,yz = zoom(psi,x,y,xv,yv,win=1)

Zoom in near (xv,yv) with a window `win=1`, which gives a 4x4 local grid.
"""
function zoom_grid(psi,x,y,xv,yv;win=1,periodic=false)
    dx,dy = Δ(x),Δ(y)
    nx,ny = size(psi)
    ix = isapprox.(x,xv,atol=dx) |> findfirst
    iy = isapprox.(y,yv,atol=dy) |> findfirst
    ixw = (ix-win):(ix+win+1)
    iyw = (iy-win):(iy+win+1)
    if !periodic
        psiz,xz,yz = psi[ixw,iyw],x[ixw],y[iyw]
    else
        ixp,iyp = mod1.(ixw,nx),mod1.(iyw,ny)
        psiz = psi[ixp,iyp] #TODO test this!
        xz = (x[ix]-win*dx):dx:(x[ix]+win*dx)
        yz = (y[iy]-win*dy):dy:(y[iy]+win*dy)
    end
    return psiz,xz,yz
end

"""
    psiz,xz,yz = zoom_interp(psi,x,y,xv,yv,win=1,nz=30)

Zoom in near (xv,yv) with a window `win=1`, which gives a 4x4 local grid.
The wavefunction psi is interpolated onto a 30x30 domain inside the local box.
"""
function zoom_interp(psi,x,y,xv,yv;win=1,nz=30,periodic=false)
    psiz, xz, yz = zoom_grid(psi,x,y,xv,yv,win=win,periodic=periodic)
    psi_itp = interpolate((xz,yz), psiz, Gridded(Linear()))
    xint,yint = LinRange(xz[1],xz[end],nz),LinRange(yz[1],yz[end],nz)
    psi_int = psi_itp(xint,yint)
    return psi_int,xint,yint
end
