function makepsi(Nv)
Lx = 300.
Ly = 150.
Nx = 1000
Ny = 500
x = collect(linspace(-Lx/2,Lx/2,Nx))
y = collect(linspace(-Ly/2,Ly/2,Ny))
dx = x[2]-x[1]
dy = y[2]-y[1]
y=y';

#randomly distributed vortices and charges
#makes sure vortices are away from edges
testvort=zeros(Nv,3)
k=1
while k<=Nv
a = -Lx/2+Lx*rand()
b = -Ly/2+Ly*rand()
σ = rand([-1,1],1)
    if (-Lx/2+dx<a<Lx/2-dx && -Ly/2+dy<b<Ly/2-dy)
        testvort[k,:] = [a b σ]
        k+=1
    end
end

testvort = sortrows(testvort)

#construct vortex wavefunction
ψ = ones(x*y)
for j=1:Nv
    ψ = makevortex(ψ,testvort[j,:],x,y)
end

return x,y,ψ,testvort
end
