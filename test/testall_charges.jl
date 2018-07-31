function testall_charges(Nv)
    x,y,psi,testvort = makepsi(Nv)
    np,nn,nt,vortices = findvortices(psi,x,y)
    vortices = remove_edgevortices(vortices,x,y)
    chargesfound = (vortices[:,3] == testvort[:,3])
    return chargesfound
  end
