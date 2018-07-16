function testallpositions(Nv)
  x,y,psi,testvort = makepsi(Nv);
  vortices = findvortices(psi,x,y);
  vortfound = checkpositions(testvort,vortices,x,y,Nv);
  return vortfound == Nv
end
