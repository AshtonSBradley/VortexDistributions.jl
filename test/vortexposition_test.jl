#find single vortex position to machine precision
using VortexDistributions, LinearAlgebra, Plots, BenchmarkTools, Revise

#make a large homogeneous wavefunction (ξ=1)
Lx = 100.; Ly = 100.
Nx = 1000; Ny = 1000
x = linspace(-Lx/2,Lx/2,Nx+1); y = linspace(-Ly/2,Ly/2,Ny+1)
x = x[1:end-1];y = y[1:end-1]
dx = x[2]-x[1]; dy = y[2]-y[1]

# Vortex at origin
x1,y1 = (0.1,0.)
vort = [x1 y1 1]

#construct vortex wavefunction (core ansatz)
psi = ones(size(x.*y')) |> complex

makeallvortices!(psi,vort,x,y,1.0)
heatmap(x,y,angle.(psi),xlabel="x",ylabel="y",transpose = true)

np,nn,nt,vortices = findvortices(psi,x,y)
vortices = remove_edgevortices(vortices,x,y)

# Random in the box interior
#place the vortices
x1,y1 = 0.8*Lx*(rand(2).-.5)
vort = [x1 y1 1]

#construct vortex wavefunction
psi = ones(size(x.*y')) |> complex

makeallvortices!(psi,vort,x,y,1.0)
heatmap(x,y,angle.(psi),xlabel="x",ylabel="y",transpose = true)

np,nn,nt,vortices = findvortices(psi,x,y)
vortices = remove_edgevortices(vortices,x,y)

xv,yv = vortices[1:2]

test1=isapprox.(x,xv,atol=dx/2)
ixv = findlast(test1)
x[ixv]
test2=isapprox.(y,yv,atol=dy/2)
iyv = findlast(test2)
y[iyv]
sum(test1)
ix = findfirst(x .== xv+dx/2)
x[ix]-dx/2
iy = findfirst(y .== yv+dy/2)
y[iy]

winhalf = 5
xwin = (ixv-winhalf):(ixv+winhalf-1)
xw = x[xwin]

ywin = (iyv-winhalf):(iyv+winhalf-1)
yw = y[ywin]

psiw = psi[xwin,ywin]

heatmap(xw,yw,angle.(psiw),xlabel="x",ylabel="y",transpose = true)



# via interpolation on window:
using Interpolations

import Interpolations:interpolate
interpolate(A::Array{Complex{Float64},2})=interpolate(real(A))+im*interpolate(imag(A))

knots = (xw,yw)
itp = interpolate(knots, psiw, Gridded(Linear()))
G(x,y) = itp(x,y)


Nf = 100

xf = LinRange(xw[winhalf],xw[winhalf+1],Nf)
yf = LinRange(yw[winhalf],yw[winhalf+1],Nf)
psif = G(xf,yf)

heatmap(xf,yf,angle.(psif),xlabel="x",ylabel="y",transpose = true)

np,nn,nt,vorticesf = findvortices(psif,xf|>Vector,yf|>Vector)
vorticesf = remove_edgevortices(vorticesf,xf,yf)

xv,yv=vorticesf[1:2]

#corezoom function (might need to disambiguiate multivortex states)
function corezoom(vortices,psi,x,y,G,winhalf=4,Nz=20)
    xv,yv = vortices[1:2]
    dx=x[2]-x[1];dy=y[2]-y[1]
    ixv = isapprox.(x,xv,atol=dx/2) |> findlast
    iyv = isapprox.(y,yv,atol=dy/2) |> findlast
    xwin = (ixv-winhalf):(ixv+winhalf-1)
    ywin = (iyv-winhalf):(iyv+winhalf-1)
    xw = x[xwin]; yw = y[ywin]
    xz = LinRange(xw[winhalf],xw[winhalf+1],Nz)
    yz = LinRange(yw[winhalf],yw[winhalf+1],Nz)
    psiz = G(xz,yz)
    np,nn,nt,vortz = findvortices(psiz,xz|>Vector,yz|>Vector)
    vortz = remove_edgevortices(vortz,xz,yz)
    return vortz,psiz,xz,yz
end

vortz,psiz,xz,yz = corezoom(vortices,psi,x,y,G)
vortz,psiz,xz,yz = corezoom(vortz,psiz,xz,yz,G)
vortz,psiz,xz,yz = corezoom(vortz,psiz,xz,yz,G)
vortz,psiz,xz,yz = corezoom(vortz,psiz,xz,yz,G)
vortz,psiz,xz,yz = corezoom(vortz,psiz,xz,yz,G)

#next iteration
dxf=xf[2]-xf[1];dyf=yf[2]-yf[1]
test1=isapprox.(xf,xv,atol=dxf/2)
ixv = findlast(test1)
xf[ixv]
test2=isapprox.(yf,yv,atol=dyf/2)
iyv = findlast(test2)
yf[iyv]
sum(test1)

winhalf = 5
xwin = (ixv-winhalf):(ixv+winhalf-1)
xw = xf[xwin]

ywin = (iyv-winhalf):(iyv+winhalf-1)
yw = yf[ywin]

xf = LinRange(xw[winhalf],xw[winhalf+1],Nf)
yf = LinRange(yw[winhalf],yw[winhalf+1],Nf)
psif = G(xf,yf)


heatmap(xf,yf,angle.(psif),xlabel="x",ylabel="y",transpose = true)

np,nn,nt,vorticesf = findvortices(psif,xf|>Vector,yf|>Vector)
vorticesf = remove_edgevortices(vorticesf,xf,yf)

xv,yv=vorticesf[1:2]
xv
x1
yv
y1
vortz[1:2]
# via second derivatives computed spectrally
# or maybe with higher order stencil derivatives?
using FFTW

angpsiw = angle.(psiw)
Lxw = xw[end]-xw[1]
Lyw = yw[end]-yw[1]

kx = (0:1:length(xw)-1)*2*pi/length(xw)|>Vector
ky = (0:1:length(xw)-1)*2*pi/length(yw)|>Vector;ky = ky'

dxpsiw = ifft(fft(angpsiw ,1).*im*(kx.*ones(length(ky))'),1)
dypsiw = ifft(fft(angpsiw ,2).*im*(ones(length(kx)).*ky),2)
dydxpsiw = ifft(fft(dxpsiw,2).*im*(ones(length(kx)).*ky),2)
dxdypsiw = ifft(fft(dypsiw,1).*im*(kx.*ones(length(ky))'),1)

vort = dxdypsiw - dydxpsiw
heatmap(xw,yw,abs.(dypsiw),xlabel="x",ylabel="y",transpose = true)
