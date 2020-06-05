using Pkg;Pkg.activate(".")
using Test, Plots, Revise, VortexDistributions
gr(transpose=true)

## field with 2 points per healing length
Nx = 200; Ny = 200
Lx = 100; Ly = Lx
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

dip = ScalarVortex([vp;vn])
periodic_dipole!(psi,dip)
heatmap(x,y,angle.(psi.ψ))

## zoom
"""
    psiz,xz,yz = zoom(psi,x,y,xc,yc,win=2)
"""
function zoom_grid(psi,x,y,xc,yc;win=2)
    dx,dy = Δ(x),Δ(y)
    ix = isapprox.(x,xc,atol=dx) |> findfirst
    iy = isapprox.(y,yc,atol=dy) |> findfirst
    ixw = (ix-win):(ix+win)
    iyw = (iy-win):(iy+win)
    xz,yz,psiz = x[ixw],y[iyw],psi[ixw,iyw]
    return psiz, xz, yz
end

function zoom_interp(psi,x,y,xc,yc;win=2,nz=30)
    psiz, xz, yz = zoom_grid(psi,x,y,xc,yc,win=win)
    psi_itp = interpolate((xz,yz), psiz, Gridded(Linear()))
    xint = LinRange(xz[1],xz[end],nz)
    yint = LinRange(yz[1],yz[end],nz)
    psi_int = psi_itp(xint,yint)
    return psi_int,xint,yint
end

## test zoom_grid
psiw,xw,yw = zoom_grid(psi.ψ,x,y,vortex_array(vp)[1:2]...)
heatmap(xw,yw,angle.(psiw))

## test zoom_interp
psi_int,xint,yint = zoom_interp(psi.ψ,x,y,vortex_array(vp)[1:2]...)
heatmap(xint,yint,angle.(psi_int))

v1 = findvortices_grid(Torus(psi_int,xint,yint),shift=true)

## test findvortices_jumps 
v2 = vort = findvortices_jumps(Torus(psi_int,xint,yint),shift=true)

psi_int,xint,yint = zoom_interp(psi_int,xint,yint,vortex_array(vp)[1:2]...)
heatmap(xint,yint,angle.(psi_int))

## definition
using Parameters, Interpolations
"""
    vortz,psiz = corezoom(vortex,ψ<:Field,winhalf=2,Nz=30)

Uses local interpolation to resolve core location.
"""
function corezoom(vortex::PointVortex,psi::T;winhalf=2,Nz=30) where T<:Field
    @unpack ψ,x,y = psi
    xv,yv,qv = vortex.x,vortex.y.vortex.q
    psiw,xw,yw = zoom_interp(psi.ψ,x,y,vortex_array(vp)[1:2]...)
    vortz = findvortices_grid(ψv,shift=true)
    vortz = remove_vortices_edge(vortz,ψv)[1]
    return vortz,ψv
end


## detect:  this must give sub-grid detection
vort = findvortices(psi)
@show vort[1]

## benchmark
using BenchmarkTools

## timing
@btime vort = findvortices(psi)
