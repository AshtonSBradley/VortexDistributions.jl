function testallpositions(Nv)
  x,y,psi,testvort=makepsi(Nv);
  vortices = findvortices(x,y,psi);
  dx = x[2]-x[1];dy=y[2]-y[1];
  vortfound = checkpositions(testvort,vortices,dx,dy,Nv);
  return vortfound == Nv
end
