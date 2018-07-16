function testallcharges(Nv)
    x,y,psi,testvort = makepsi(Nv);
    vortices = findvortices(psi,x,y);
    vortfound = checkpositions(testvort,vortices,x,y,Nv);
    chargesfound = (vortices[:,3] == testvort[:,3])
    return chargesfound
  end
