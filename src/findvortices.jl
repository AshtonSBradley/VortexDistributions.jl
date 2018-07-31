"""
   vortices = findvortices(ψ,x,y)

Locates vortices as 2π phase windings around plaquettes on a cartesian spatial field.

Requires a 2D wavefunction ψ(x,y) specified on a cartesian grid.
"""

function findvortices(ψ,x,y)
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

   ixp,iyp,vp = findnz(diffx.>0.)
   xp = x[ixp]; yp = y[iyp]
   np = length(vp)

   ixn,iyn,vn = findnz(diffx.<0.)
   xn = x[ixn]; yn = y[iyn];
   nn = length(vn)
   nt = np + nn

   vortices = [xn yn -vn; xp yp vp] |> sortrows
   return np,nn,nt,vortices
end
