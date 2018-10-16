"""
`vortz,psiz,xz,yz = zcorezoom(vortices,psi,x,y,winhalf=2,Nz=30)`

Uses local interpolation to resolve core location to ~ 5 figures"""
function corezoom(vortices,psi,x,y,winhalf=2,Nz=30)
    println("this one")
    xv,yv = vortices[1:2]
    dx=x[2]-x[1];dy=y[2]-y[1]
    ixv = isapprox.(x,xv,atol=dx) |> findlast
    iyv = isapprox.(y,yv,atol=dy) |> findlast
    ixwin = (ixv-winhalf):(ixv+winhalf-1)
    iywin = (iyv-winhalf):(iyv+winhalf-1)
    xw = x[ixwin]; yw = y[iywin]; psiw = psi[ixwin,iywin]
    xz = LinRange(xw[1],xw[end],Nz)
    yz = LinRange(yw[1],yw[end],Nz)
    knots = (xw,yw)
    itp = interpolate(knots, psiw, Gridded(Linear()))
    psiz = itp(xz,yz)
    np,nn,nt,vortz = findvortices(psiz,xz|>Vector,yz|>Vector)
    nt,np,nn,vortz = remove_edgevortices(vortz,xz|>Vector,yz|>Vector)
    return vortz,psiz,xz,yz
end
