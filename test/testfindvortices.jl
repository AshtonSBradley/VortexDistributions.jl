using Plots, LinearAlgebra, Revise, VortexDistributions

# test on randomly placed vortices
Lx = 300.; Ly = 150.
Nx = 1000; Ny = 500
x = linspace(-Lx/2,Lx/2,Nx)
y = linspace(-Ly/2,Ly/2,Ny)
ξ = 1.0

Nv = 10
testvort = randomvortices(x,y,Nv)

ψ = ones(size(x.*y')) |> complex
makeallvortices!(ψ,testvort,x,y,ξ)

ϕ = angle.(ψ)
heatmap(x,y,ϕ,transpose=true)
xlabel!("x");ylabel!("y")

nt,np,nn,vortices = findvortices(ψ,x,y)
nt,np,nn,vortices = remove_edgevortices(vortices,x,y)
chargesfound = (vortices[:,3] == testvort[:,3])

vortfound = checkvortexlocations(testvort,vortices,x,y,Nv)

vortfound == Nv

vortices



testvort
# type version
psiv = Torus(x,y,ψ)

nt,np,nn,vortices = findvortices(psiv)
nt,np,nn,vortices = removeedgevortices(vortices,x,y)

chargesfound = (vortices[:,3] == testvort[:,3])

vortfound = checkvortexlocations(testvort,vortices,x,y,Nv)

vortfound == Nv


# test on saved data
using Plots, JLD2, Parameters, Interpolations

@load "./examples/one_frame.jld2"

gr(transpose=false)
heatmap(abs2.(ψ1'))
heatmap(x,y,angle.(ψ1'))

xind = 100:400 #100:300 for no vortices
yind = 200:800
xwin = x[xind]
ywin = y[yind]
ψwin = ψ1[xind,yind]
heatmap(xwin,ywin,abs2.(ψwin'))
heatmap(xwin,ywin,angle.(ψwin'))

# make some synthetic data
psi2 = randn(size(ψ1))+im*randn(size(ψ1))
heatmap(x,y,angle.(psi2'))
ψwin = psi2[xind,yind]

psi = Torus(x,y,ψ1)
psiw = Torus(xwin,ywin,ψwin)

nt,np,nn,vortices = findvorticesinterp(psiw)

nt,np,nn,vortices = findvorticesinterp(psiw)

function makepsi(Nv,Lx=200,Ly=200,Nx=400,Ny=400)

x = linspace(-Lx/2,Lx/2,Nx)
y = linspace(-Ly/2,Ly/2,Ny)

#randomly distributed vortices and charges
#make sure vortices are away from edges
#make sure vortices don't have same x-coordinates

testvort = randomvortices(x,y,Nv)

#construct vortex wavefunction
ψ = ones(size(x.*y')) |> complex
makeallvortices!(ψ,testvort,x,y,.1)

return x,y,ψ,testvort
end

x,y,ψ3,testvort = makepsi(2)

makeallvortices!(ψ3,testvort,x,y)

psi3 = Torus(x,y,ψ3)
nt,np,nn,vortices = findvortices(psi3)
nt,np,nn,vortices = findvorticesgrid(psi3)
