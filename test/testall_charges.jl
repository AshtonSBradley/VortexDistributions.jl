function testall_charges(Nv)
    x,y,psi,testvort = makepsi(Nv);
    vortices = findvortices(psi,x,y);
    chargesfound = (vortices[:,3] == testvort[:,3])
    return chargesfound
  end
