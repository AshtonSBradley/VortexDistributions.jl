function findwhere(A)
    I = findall(!iszero,A)
    v = A[I]
    ix = [I[i][1] for i in eachindex(I)]
    iy = [I[i][2] for i in eachindex(I)]
    return ix,iy,v
end

"""
    vortices = removeedgevortices(vort::Array{PointVortex,1},x,y,edge=1)

Strip artifact edgevortices arising from periodic phase differencing.
"""
function removeedgevortices(vort::Array{PointVortex,1},psi::Field,edge=1)
    @unpack x,y = psi; dx=x[2]-x[1]; dy=y[2]-y[1]
    keep = []
    for j = 1:length(vort)
        xi,yi,qi = rawData(vort[j])
        xedge = isapprox(xi,x[1],atol=edge*dx) || isapprox(xi,x[end],atol=edge*dx)
        yedge = isapprox(yi,y[1],atol=edge*dy) || isapprox(yi,y[end],atol=edge*dy)
        not_edge = !(xedge || yedge)
        not_edge && push!(keep,j)
    end
    return vort[keep]
end
