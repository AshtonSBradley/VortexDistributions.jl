using Pkg;Pkg.activate(".")
using Test, Plots, BenchmarkTools, Revise
using VortexDistributions
gr(transpose=true)

## field with 2 points per healing length
Nx = 400; Ny = 400
Lx = 200; Ly = Lx
x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, Ny+1)[1:end-1];
psi0 = one.(x*y') |> complex;

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

## import for debugging
using Parameters, Interpolations
import VortexDistributions:findvortices, zoom_interp, zoom_grid

function findvortices(psi::Field;periodic=false)
    @unpack ψ,x,y = psi
    vort = findvortices_grid(psi)
    if !periodic
        vort = remove_vortices_edge(vort,psi) #not periodic: remove
    end

    for (j,vortex) in enumerate(vort)
        v = try
        psi_int,xint,yint = zoom_interp(ψ,x,y,vortex.xv,vortex.yv,periodic=periodic) #periodic: peridic indices here
        v1 = findvortices_grid(Torus(psi_int,xint,yint))
        vint = remove_vortices_edge(v1,Torus(psi_int,xint,yint))[1]
        psi_int,xint,yint = zoom_interp(psi_int,xint,yint,vint.xv,vint.yv)
        v1 = findvortices_grid(Torus(psi_int,xint,yint))
        vint = remove_vortices_edge(v1,Torus(psi_int,xint,yint))[1]
        catch nothing
        end
        v != nothing && (vort[j] = v)    # NOTE fallback to grid if zoom fails
    end
    if periodic
        Lx,Ly = last(x)-first(x),last(y)-first(y)
        vdat = vortex_array(vort)
        xx,yy,qq = vdat[:,1],vdat[:,2],vdat[:,3]
        @. xx[xx>Lx/2] -= Lx
        @. xx[xx<-Lx/2] += Lx
        @. yy[yy>Ly/2] -= Ly
        @. yy[yy<-Lx/2] += Ly
        vort = PointVortex.(xx,yy,qq)
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

## test !periodic detection
vort = findvortices(psi)
@test length(vort) == 0

# @btime findvortices(psi)

## test periodic detection
vort = findvortices(psi,periodic=true)
@test length(vort) == 2

# @btime findvortices(psi,periodic=true)

## single vortex test
psi = Torus(copy(psi0),x,y)

xp = 1.133
yp = 2.787
vp = PointVortex(xp,yp,1)

sp = ScalarVortex(vp)

vortex!(psi,sp)
vfound = findvortices(psi);@show vortex_array(vfound)

## test zoom window requirements and stability for single vortex
# round one:

vort = findvortices_grid(psi)
vort = remove_vortices_edge(vort,psi)

psi_int,xint,yint = zoom_interp(psi.ψ,psi.x,psi.y,vort[1].xv,vort[1].yv,periodic=true)

v1 = findvortices_grid(Torus(psi_int,xint,yint))
vint = remove_vortices_edge(v1,Torus(psi_int,xint,yint))[1]

@show vint
heatmap(xint,yint,angle.(psi_int))
scatter!([vint.xv],[vint.yv],label=false)

## round two
psi_int,xint,yint = zoom_interp(psi_int,xint,yint,vint.xv,vint.yv)
v2 = findvortices_grid(Torus(psi_int,xint,yint))
vint = remove_vortices_edge(v2,Torus(psi_int,xint,yint))[1]
@show vint
heatmap(xint,yint,angle.(psi_int))
scatter!([vint.xv],[vint.yv],label=false)

## test zoom stages on failing cases in real data
using JLD2
loadpath="/Users/abradley/Dropbox/Julia/Vortices - simple results"
f1 = joinpath(loadpath,"nv20-t_0.jld2")
f2 = joinpath(loadpath,"nv20-t_200.jld2")



## load good data
@load f2 𝛹200 x y
psi200 = Torus(𝛹200,x,y)
vfound = findvortices(psi200,periodic=true)
p1=heatmap(x,y,angle.(psi200.ψ))
vdata = vortex_array(vfound)



for j in 1:size(vdata)[1]
    scatter!(p1,[vdata[j,1]],[vdata[j,2]],label=false,color=(vdata[j,3]==1 ? :green : :blue),ms=4,alpha=.6,markerstrokecolor=(vdata[j,3]==1 ? :green : :blue))
end
p1

## load offending data
@load f1 𝛹0 x y
psi = Torus(𝛹0,x,y)
p1=heatmap(x,y,angle.(psi.ψ))
vfound = findvortices(psi,periodic=true)
vdata = vortex_array(vfound)

for j in 1:size(vdata)[1]
    scatter!(p1,[vdata[j,1]],[vdata[j,2]],label=false,color=(vdata[j,3]==1 ? :green : :blue),ms=4,alpha=.6,markerstrokecolor=(vdata[j,3]==1 ? :green : :blue))
end
p1

## choose vortex that seem
