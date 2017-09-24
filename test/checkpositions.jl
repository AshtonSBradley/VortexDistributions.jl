function checkpositions(testvort,vortices,dx,dy,Nv)
#check detection to grid resolution
vortfound = 0
for j=1:Nv
    if abs.(testvort[j,1]-vortices[j,1])<dx && abs.(testvort[j,2]-vortices[j,2])<dy
        vortfound+=1
    end
end
return vortfound
end
