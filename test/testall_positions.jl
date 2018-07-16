function testall_positions(Nv)
    x,y,psi,testvort = makepsi(Nv)
    vortices = findvortices(psi,x,y)
    vortices = remove_edgevortices(vortices,x,y)
    vortfound = checkvortexlocations(testvort,vortices,x,y,Nv)
    return vortfound == Nv
end
