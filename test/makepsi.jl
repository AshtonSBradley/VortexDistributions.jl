function makepsi(Nv,Lx,Ly,Nx,Ny)

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
