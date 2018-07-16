function testall_charges(Nv)
    x,y,psi,testvort = makepsi(Nv);
    psi = edgemask(psi,x,y)
    vortices = findvortices(psi,x,y);
    chargesfound = (vortices[:,3] == testvort[:,3])
    return chargesfound
  end
