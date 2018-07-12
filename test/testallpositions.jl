function testallpositions(Nv)
  x,y,psi,testvort = makepsi(Nv);
  vortices = findvortices(x,y,psi);
  vortfound = checkpositions(testvort,vortices,x,y,Nv);
  return vortfound == Nv
end
