Lx = 300.
Ly = 150.
Nx = 1000
Ny = 500
x = linspace(-Lx/2,Lx/2,Nx) |> collect
y = linspace(-Ly/2,Ly/2,Ny) |> collect
dx = x[2]-x[1]
dy = y[2]-y[1]

#randomly distributed vortices and charges
#makes sure vortices are away from edges
Nv = 5
testvort=zeros(Nv,3)
k=1
while k<=Nv
a = -Lx/2 + Lx*rand()
b = -Ly/2 + Ly*rand()
σ = rand([-1,1],1)
    if (-Lx/2 + dx < a < Lx/2 - dx && -Ly/2 + dy < b < Ly/2 - dy)
        testvort[k,:] = [a b σ]
        global k+=1
    end
end

testvort = sortslices(testvort,dims=1)

#construct vortex wavefunction
ψ = ones(size(x.*y')) |> complex
for j=1:Nv
    makevortex!(ψ,testvort[j,:],x,y,.1)
end
