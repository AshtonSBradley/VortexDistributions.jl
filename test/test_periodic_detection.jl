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
