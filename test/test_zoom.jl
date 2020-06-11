using Pkg;Pkg.activate(".")
using Test, Plots, BenchmarkTools, Revise, VortexDistributions
gr(transpose=true)

## zoom
using Interpolations, Parameters
"""
    psiz,xz,yz = zoom(psi,x,y,xv,yv,win=1)

Zoom in near (xv,yv) with a window `win=1`, which gives a 4x4 local grid.
"""
function zoom_grid(psi,x,y,xv,yv;win=1)
    dx,dy = Δ(x),Δ(y)
    ix = isapprox.(x,xv,atol=dx) |> findfirst
    iy = isapprox.(y,yv,atol=dy) |> findfirst
    ixw = (ix-win):(ix+win+1)
    iyw = (iy-win):(iy+win+1)
    psiz,xz,yz = psi[ixw,iyw],x[ixw],y[iyw]
    return psiz,xz,yz
end

"""
    psiz,xz,yz = zoom_interp(psi,x,y,xv,yv,win=1,nz=30)

Zoom in near (xv,yv) with a window `win=1`, which gives a 4x4 local grid.
The wavefunction psi is interpolated onto a 30x30 domain inside the local box.
"""
function zoom_interp(psi,x,y,xv,yv;win=1,nz=30)
    psiz, xz, yz = zoom_grid(psi,x,y,xv,yv,win=win)
    psi_itp = interpolate((xz,yz), psiz, Gridded(Linear()))
    xint,yint = LinRange(xz[1],xz[end],nz),LinRange(yz[1],yz[end],nz)
    psi_int = psi_itp(xint,yint)
    return psi_int,xint,yint
end

## field with 2 points per healing length
Nx = 400; Ny = 400
Lx = 200; Ly = Lx
x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, Ny+1)[1:end-1];
psi0 = one.(x*y') |> complex

## make vortex dipole
# Periodic boundary conditions
psi = Torus(copy(psi0),x,y)

xp = 1.7
yp = 2.1
vp = PointVortex(xp,yp,1)

xn = 3.25
yn = 10.
vn = PointVortex(xn,yn,-1)

dip = ScalarVortex([vp;vn])
periodic_dipole!(psi,dip)

heatmap(x,y,angle.(psi.ψ))

vfound = findvortices(psi)

## single vortex
# Periodic boundary conditions
psi = Torus(copy(psi0),x,y)

xp = 1.13
yp = 2.78
vp = PointVortex(xp,yp,1)

sp = ScalarVortex(vp)

vortex!(psi,sp)

heatmap(x,y,angle.(psi.ψ))
## test zoom_grid

psiw,xw,yw = zoom_grid(psi.ψ,x,y,vortex_array(vp)[1:2]...)
heatmap(xw,yw,angle.(psiw))

## test zoom_interp
# TODO mystery of why shift has no effect but to make x grid resolved...
# make sure zoom stability is exact
# test findvortices_jumps has periodic correction that will need to be rethought.


@info "start zoom for known (?) vortex"
psi_int,xint,yint = zoom_interp(psi.ψ,x,y,vortex_array(vp)[1:2]...)
v1 = findvortices_grid(Torus(psi_int,xint,yint))
vint = remove_vortices_edge(v1,Torus(psi_int,xint,yint))[1]
@show vint

# second iteration

psi_int,xint,yint = zoom_interp(psi_int,xint,yint,vint.xv,vint.yv)
v2 = findvortices_grid(Torus(psi_int,xint,yint))
vint = remove_vortices_edge(v2,Torus(psi_int,xint,yint))[1]
@show vint
heatmap(xint,yint,angle.(psi_int))
scatter!([vint.xv],[vint.yv],label=false)

# conclusion seems to be that 2 iterations is giving high accuracy and stability

## test findvortices

@btime vfound = findvortices(psi)
