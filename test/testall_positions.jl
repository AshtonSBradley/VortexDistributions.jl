function testall_positions(Nv)
  x,y,psi,testvort = makepsi(Nv);
  psi = edgemask(psi,x,y)
  vortices = findvortices(psi,x,y);
  vortfound = checkvortexlocations(testvort,vortices,x,y,Nv);
  return vortfound == Nv
end
