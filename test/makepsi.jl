function makepsi(Nv)
Lx = 300.
Ly = 150.
Nx = 1000
Ny = 500
x = linspace(-Lx/2,Lx/2,Nx)
y = linspace(-Ly/2,Ly/2,Ny)

#randomly distributed vortices and charges
#make sure vortices are away from edges
#make sure vortices don't have same x-coordinates

testvort = randomvortices(x,y,Nv)

#construct vortex wavefunction
ψ = ones(size(x.*y')) |> complex
for j=1:Nv
    makevortex!(ψ,testvort[j,:],x,y,.1)
end

return x,y,ψ,testvort
end
