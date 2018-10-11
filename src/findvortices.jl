"""
   `np,nn,nt,vortices = findvortices(ψ,x,y)`

Locates vortices as 2π phase windings around plaquettes on a cartesian spatial grid. Uses an optimized plaquette method.

Requires a 2D wavefunction ψ(x,y) on a cartesian grid specified by vectors x, y.

`np` - number of positive vortices

`nn` - number of negative vortices

`nt` - total number of vortices

`vortices` - array of vortex coordinates `xv,yv` and circulations `cv`. Each row is of the form `[xv, yv, cv]`, and the array is sorted into lexical order according to the `xv` coordinates
"""

function findvortices(ψ,x,y;geometry="torus")
@assert typeof(x)==Array{Float64,1}
@assert typeof(y)==Array{Float64,1}
@assert typeof(ψ)==Array{Complex{Float64},2}

   phase = angle.(ψ)
   diffx = countphasejumps(phase,1)
   diffy = countphasejumps(phase,2)
   #circ = diffx .- circshift(diffx,(0,1)) .- diffy .+ circshift(diffy,(1,0))
   #diffx .-= circshift(diffx,(0,1))
   #diffx .-= diffy
   #diffx .+= circshift(diffy,(1,0))

   #use in-place memory recycling
   circshift!(phase,diffx,(0,1))
   diffx .-= phase
   diffx .-= diffy
   circshift!(phase,diffy,(1,0))
   diffx .+= phase

   ixp,iyp,vp = findnonzero(diffx.>0.)
   xp = x[ixp]; yp = y[iyp]
   np = length(vp)

   ixn,iyn,vn = findnonzero(diffx.<0.)
   xn = x[ixn]; yn = y[iyn];
   nn = length(vn)

   nt = np + nn

   #offest for plaquette size (TODO: add accurate solver a-la Ricardo)
   dx = x[2]-x[1]; dy = y[2] - y[1]
   xp += -dx/2; yp += -dx/2; xn += -dx/2; yn += -dx/2
   
if geometry == "sphere"
    #find polar winding
    vortices = [xn yn -vn; xp yp vp]

    # to do: optimize this step at start of method to avoid extra angle call:
    windvals = countphasejumps(angle.(ψ),2)

    #test north pole for winding. - sign means polar vortex co-rotating with w > 0
    w1 = sum(windvals[1,:])
    (sign(w1) != 0) && (vortices = [vortices; 0.0 0.0 -w1])

    #test south pole for winding. + sign means co-rotating with w > 0
    w2 = sum(windvals[end,:])
    (sign(w2) != 0) && (vortices = [vortices; pi 0.0 w2])

elseif geometry == "torus"
    vortices = [xn yn -vn; xp yp vp]
end
    vortices = sortslices(vortices,dims=1)
return nt,np,nn,vortices
end
