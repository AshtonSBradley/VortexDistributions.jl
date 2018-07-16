function remove_edgevortices(vortices,x,y)
#eliminate edge vortices (makepsi creates non-periodic wavefunction)
dx=x[2]-x[1];dy=y[2]-y[1];
Nfound,_ = size(vortices)
keep = []
for j = 1:Nfound
    xi = vortices[j,1]; yi = vortices[j,2];
    xedge = isapprox(xi,x[1],atol=dx) || isapprox(xi,x[end],atol=dx)
    yedge = isapprox(yi,y[1],atol=dy) || isapprox(yi,y[end],atol=dy)
    not_edge = !(xedge || yedge)
    not_edge && push!(keep,j)
end
vortices=vortices[keep,:]
    return vortices
end
