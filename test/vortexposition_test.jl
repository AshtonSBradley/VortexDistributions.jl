#find single vortex position to machine precision
using Pkg
pkg"activate ."

using VortexDistributions, LinearAlgebra, Plots, BenchmarkTools, Revise

# via interpolation on window:
# using Interpolations
# Interpolations.interpolate(A::Array{Complex{Float64},2})=interpolate(real(A))+im*interpolate(imag(A))


# make a large homogeneous wavefunction (ξ=1)
Lx = 100.; Ly = 150.
Nx = 300; Ny = 400
x = linspace(-Lx/2,Lx/2,Nx+1); y = linspace(-Ly/2,Ly/2,Ny+1)
x = x[1:end-1];y = y[1:end-1]
dx = x[2]-x[1]; dy = y[2]-y[1]

# Vortex at origin
x1,y1 = (0.,0.)
vort = [x1 y1 1]

#construct vortex wavefunction (core ansatz)
psi = ones(size(x.*y')) |> complex
makeallvortices!(psi,vort,x,y,1.0)
heatmap(x,y,angle.(psi),xlabel="x",ylabel="y",transpose = true)

nt,np,nn,vortices = findvortices(psi,x,y,false)
nt,np,nn,vortices = remove_edgevortices(vortices,x,y)


# Random in the box interior
#place the vortices
x1,y1 = Lx*(rand(2).-.5)
vort = [x1 y1 1]
vort[1:2]
#construct vortex wavefunction
psi = ones(size(x.*y')) |> complex
makeallvortices!(psi,vort,x,y,1.0)
heatmap(x,y,angle.(psi),xlabel="x",ylabel="y",transpose = true)

nt,np,nn,vortices = findvortices(psi,x,y,false)
nt,np,nn,vortices = remove_edgevortices(vortices,x,y)

xv,yv = vortices[1:2]



# ixv = isapprox.(x,xv,atol=dx/2) |> findlast
# iyv = isapprox.(y,yv,atol=dy/2) |> findlast

# winhalf = 6
# xwin = (ixv-winhalf):(ixv+winhalf-1)
# xw = x[xwin]
#
# ywin = (iyv-winhalf):(iyv+winhalf-1)
# yw = y[ywin]
#
# psiw = psi[xwin,ywin]
#
# heatmap(xw,yw,angle.(psiw),xlabel="x",ylabel="y",transpose = true)

#knots = (xw,yw)
#itp = interpolate(knots, psiw, Gridded(Linear()))

# THIS BLOCK
# knots = (x,y)
# itp = interpolate(knots, psi, Gridded(Linear()))
#
# corezoom function (might need to disambiguiate multivortex states)
# function corezoom(vortices,psi,x,y,winhalf=2,Nz=30)
#     xv,yv = vortices[1:2]
#     dx=x[2]-x[1];dy=y[2]-y[1]
#     ixv = isapprox.(x,xv,atol=dx) |> findlast
#     iyv = isapprox.(y,yv,atol=dy) |> findlast
#     ixwin = (ixv-winhalf):(ixv+winhalf-1)
#     iywin = (iyv-winhalf):(iyv+winhalf-1)
#     xw = x[ixwin]; yw = y[iywin]; psiw = psi[ixwin,iywin]
#     xz = LinRange(xw[1],xw[end],Nz)
#     yz = LinRange(yw[1],yw[end],Nz)
#     knots = (xw,yw)
#     itp = interpolate(knots, psiw, Gridded(Linear()))
#     psiz = itp(xz,yz)
#     np,nn,nt,vortz = findvortices(psiz,xz|>Vector,yz|>Vector)
#     nt,np,nn,vortz = remove_edgevortices(vortz,xz|>Vector,yz|>Vector)
#     return vortz,psiz,xz,yz
# end
winhalf = 2
Nz = 30
xv,yv=vortices[1:2]
dx=x[2]-x[1];dy=y[2]-y[1]
ixv = isapprox.(x,xv,atol=dx) |> findlast
iyv = isapprox.(y,yv,atol=dy) |> findlast
ixwin = (ixv-winhalf):(ixv+winhalf-1)
iywin = (iyv-winhalf):(iyv+winhalf-1)
xw = x[ixwin];yw = y[iywin]; psiw = psi[ixwin,iywin]
xz = LinRange(xw[1],xw[end],Nz)
yz = LinRange(yw[1],yw[end],Nz)
knots = (xw,yw)
itp = interpolate(knots,psiw,Gridded(Linear()))
psiz = itp(xz,yz)
nt,np,nn,vortz = findvortices_grid(psiz,xz|>Vector,yz|>Vector)
nt,np,nn,vortz = remove_edgevortices(vortz,xz|>Vector,yz|>Vector)





vortz,psiz,xz,yz = corezoom(vortices,psi,x,y)
vortz,psiz,xz,yz = corezoom(vortz,psiz,xz,yz)
vortz,psiz,xz,yz = corezoom(vortz,psiz,xz,yz)
#heatmap(xz,yz,angle.(psiz),xlabel="x",ylabel="y",transpose = true)

xv,yv = vortz[1:2]

show(xv)
show(x1)

show(yv)
show(y1)
#WORKS!

#DOESN'T:
#testing in package
using Pkg
pkg"activate ."

using VortexDistributions
include("makepsi.jl")
Lx = 100.
Ly = 100.
Nx = 300
Ny = 300
Nv = 2
x,y,psi,testvortices = makepsi(Nv,Lx,Ly,Nx,Ny)

nt,np,nn,vortices = findvortices(psi,x,y)
nt,np,nn,vortices = remove_edgevortices(vortices,x,y)




# for j in 1:nt
#     vortex = vortices[j,:]
#     vortz,psiz,xz,yz = corezoom(vortex,psi,x,y)
#     vortz,psiz,xz,yz = corezoom(vortz,psiz,xz,yz)
#     vortz,psiz,xz,yz = corezoom(vortz,psiz,xz,yz)
#     vortices[j,1:2] = vortz[1:2]
# end
