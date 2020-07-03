"""
`vortices = remove_edgevortices(vortices,x,y,edge=1)`

Strips edgevortices due to periodic phase differencing."""
function remove_vortices_edge(vortices,x,y,edge=1)
    dx=x[2]-x[1];dy=y[2]-y[1];
    Nfound,_ = size(vortices)
    keep = []
    for j = 1:Nfound
        xi = vortices[j,1]; yi = vortices[j,2];
        xedge = isapprox(xi,x[1],atol=edge*dx) || isapprox(xi,x[end],atol=edge*dx)
        yedge = isapprox(yi,y[1],atol=edge*dy) || isapprox(yi,y[end],atol=edge*dy)
        not_edge = !(xedge || yedge)
        not_edge && push!(keep,j)
    end
    vortices = vortices[keep,:]
    nt = size(vortices)[1]
    np = sum(vortices[:,3].>0)
    nn = nt - np
    return nt,np,nn,vortices
end
