#find single vortex position to machine precision
using VortexDistributions, LinearAlgebra, Plots, BenchmarkTools, Revise

# via interpolation on window:
using Interpolations
Interpolations.interpolate(A::Array{Complex{Float64},2})=interpolate(real(A))+im*interpolate(imag(A))


#make a large homogeneous wavefunction (ξ=1)
Lx = 100.; Ly = 100.
Nx = 1000; Ny = 1000
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

np,nn,nt,vortices = findvortices(psi,x,y)
vortices = remove_edgevortices(vortices,x,y)

# Random in the box interior
#place the vortices
x1,y1 = Lx*(rand(2).-.5)
vort = [x1 y1 1]

#construct vortex wavefunction
psi = ones(size(x.*y')) |> complex

makeallvortices!(psi,vort,x,y,1.0)
heatmap(x,y,angle.(psi),xlabel="x",ylabel="y",transpose = true)

np,nn,nt,vortices = findvortices(psi,x,y)
vortices = remove_edgevortices(vortices,x,y)

xv,yv = vortices[1:2]
ixv = isapprox.(x,xv,atol=dx/2) |> findlast
iyv = isapprox.(y,yv,atol=dy/2) |> findlast

winhalf = 4
xwin = (ixv-winhalf):(ixv+winhalf-1)
xw = x[xwin]

ywin = (iyv-winhalf):(iyv+winhalf-1)
yw = y[ywin]

psiw = psi[xwin,ywin]

heatmap(xw,yw,angle.(psiw),xlabel="x",ylabel="y",transpose = true)

knots = (xw,yw)
itp = interpolate(knots, psiw, Gridded(Linear()))

#corezoom function (might need to disambiguiate multivortex states)
function corezoom(vortices,psi,x,y,itp,knots,winhalf=4,Nz=20)
    xv,yv = vortices[1:2]
    dx=x[2]-x[1];dy=y[2]-y[1]
    ixv = isapprox.(x,xv,atol=dx/2) |> findlast
    iyv = isapprox.(y,yv,atol=dy/2) |> findlast
    ixwin = (ixv-winhalf):(ixv+winhalf-1)
    iywin = (iyv-winhalf):(iyv+winhalf-1)
    xw = x[ixwin]; yw = y[iywin]; psiw = psi[ixwin,iywin]
    xz = LinRange(xw[1],xw[end],Nz)
    yz = LinRange(yw[1],yw[end],Nz)
    #knots = (xw,yw)
    #itp = interpolate(knots, psi, Gridded(Linear()))
    psiz = itp(xz,yz)
    np,nn,nt,vortz = findvortices(psiz,xz|>Vector,yz|>Vector)
    vortz = remove_edgevortices(vortz,xz|>Vector,yz|>Vector)
    return vortz,psiz,xz,yz
end


winhalf = 4
Nz = 50
vortz,psiz,xz,yz = corezoom(vortices,psi,x,y,itp,knots,winhalf,Nz)
heatmap(xz,yz,angle.(psiz),xlabel="x",ylabel="y",transpose = true)
vortz,psiz,xz,yz = corezoom(vortz,psiz,xz,yz,itp,knots,winhalf,Nz)
heatmap(xz,yz,angle.(psiz),xlabel="x",ylabel="y",transpose = true)
vortz,psiz,xz,yz = corezoom(vortz,psiz,xz,yz,itp,knots,winhalf,Nz)
heatmap(xz,yz,angle.(psiz),xlabel="x",ylabel="y",transpose = true)

xv,yv = vortz[1:2]

xv
x1

yv
y1
