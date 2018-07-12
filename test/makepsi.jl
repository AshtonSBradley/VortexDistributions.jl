function makepsi(Nv)
Lx = 300.
Ly = 150.
Nx = 1000
Ny = 500
x = linspace(-Lx/2,Lx/2,Nx) |> collect
y = linspace(-Ly/2,Ly/2,Ny) |> collect
dx = x[2]-x[1]; dy = y[2]-y[1]

#randomly distributed vortices and charges
#make sure vortices are away from edges
#make sure vortices don't have same x-coordinates

testvort = zeros(Nv,3)
k = 1
while k <= Nv
a = rand(x[2:end-1])
b = rand(y[2:end-1])
σ = rand([-1,1],1)

        xdist = abs.(a-testvort[:,1])
            if minimum(xdist) > 2*dx
                testvort[k,:] = [a b σ]
                k += 1
            end

end

testvort = sortrows(testvort)

#construct vortex wavefunction
ψ = complex(ones(x.*y'))
for j=1:Nv
    makevortex!(ψ,testvort[j,:],x,y,.1)
end

return x,y,ψ,testvort
end
