using VortexDistributions, PyPlot, LinearAlgebra, Revise

# create some grids
Nv = 50
Lx = 300.; Ly = 150.
Nx = 1000; Ny = 500
x = linspace(-Lx/2,Lx/2,Nx)
y = linspace(-Ly/2,Ly/2,Ny)

testvort = randomvortices(x,y,Nv)

ξ = 1.0

ψ = ones(size(x.*y')) |> complex

#=
for j=1:Nv
    makeallvortex!(ψ,testvort[j,:],x,y,ξ)
end
=#
makeallvortices!(ψ,testvort,x,y,ξ)

psi = zeros(Ny,Nx) |> complex
transpose!(psi,ψ)
#psi = reverse(psi,dims=1)
ϕ = angle.(psi)
pcolormesh(x,y,ϕ)
xlabel("x");ylabel("y")
colorbar()

nt,np,nn,vortices = findvortices(ψ,x,y)
vortices = remove_edgevortices(vortices,x,y)
