function checkvortexlocations(testvort,vortices,x,y,Nv)
#check detection to grid resolution
dx = x[2] - x[1]; dy = y[2] - y[1]
vortfound = 0
for j=1:Nv
    if (isapprox(testvort[j,1],vortices[j,1],atol=dx) && isapprox(testvort[j,2],vortices[j,2],atol=dy))
        vortfound+=1
    end
end
return vortfound
end
