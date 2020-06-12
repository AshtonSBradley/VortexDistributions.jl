using Pkg;Pkg.activate(".")
using Test, Plots, BenchmarkTools, Revise
using VortexDistributions
gr(transpose=true)

## field with 2 points per healing length
Nx = 400; Ny = 400
Lx = 200; Ly = Lx
x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, Ny+1)[1:end-1];
psi0 = one.(x*y') |> complex

## make a field with a dipole
# Periodic boundary conditions
psi = Torus(copy(psi0),x,y)

xp = 99.5
yp = 0
xn = xp
yn = -50
vp = PointVortex(xp,yp,1)
vn = PointVortex(xp,yn,-1)
dip = ScalarVortex([vp;vn])

periodic_dipole!(psi,dip)

heatmap(x,y,angle.(psi.ψ))

## periodic
using Parameters, Interpolations
import VortexDistributions:findvortices, zoom_interp, zoom_grid
function findvortices(psi::Field)
    @unpack ψ,x,y = psi
    vort = findvortices_grid(psi)
    # vort = remove_vortices_edge(vort,psi) #periodic: don't remove

    for (j,vortex) in enumerate(vort)
        v = try
        psi_int,xint,yint = zoom_interp(ψ,x,y,vortex.xv,vortex.yv,periodic=true) #periodic: peridic indices here
        v1 = findvortices_grid(Torus(psi_int,xint,yint))
        vint = remove_vortices_edge(v1,Torus(psi_int,xint,yint))[1]
        psi_int,xint,yint = zoom_interp(psi_int,xint,yint,vint.xv,vint.yv)
        v1 = findvortices_grid(Torus(psi_int,xint,yint))
        vint = remove_vortices_edge(v1,Torus(psi_int,xint,yint))[1]
        catch nothing
        end
        vort[j] = v     # NOTE fallback to grid if zoom fails
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
        ixp,iyp = mod1.(ixw,nx),mod1.(iyw,ny) # periodic indices
        psiz = psi[ixp,iyp]
        xz = (x[ix]-win*dx):dx:(x[ix]+(win+1)*dx) # local x vector
        yz = (y[iy]-win*dy):dy:(y[iy]+(win+1)*dy) # local y vector
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

# @btime vort = findvortices(psi)

vort = findvortices(psi)

## detect
vort = findvortices_grid(psi)

psi_int,xint,yint = zoom_interp(psi.ψ,psi.x,psi.y,vort[1].xv,vort[1].yv,periodic=true)

heatmap(xint,yint,angle.(psi_int))

v1 = findvortices_grid(Torus(psi_int,xint,yint))
vint = remove_vortices_edge(v1,Torus(psi_int,xint,yint))[1]

## round two
psi_int,xint,yint = zoom_interp(psi_int,xint,yint,vint.xv,vint.yv)
v2 = findvortices_grid(Torus(psi_int,xint,yint))
vint = remove_vortices_edge(v2,Torus(psi_int,xint,yint))[1]
@show vint
heatmap(xint,yint,angle.(psi_int))
