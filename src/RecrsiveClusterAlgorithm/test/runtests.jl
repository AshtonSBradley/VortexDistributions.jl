using Pkg
Pkg.activate(".")
using Test, LinearAlgebra, Plots, Revise, RecursiveClusterAlgorithm
gr(grid=false,legend=false)

## get_dipoles
# some vortices
xi = [-0.549 1.160 0.335 -1.773]
yi = [0.123 -0.603 -0.079 0.220]

scatter(xi[1:2],yi[1:2],m=:circle,c=:red)
scatter!(xi[3:4],yi[3:4],m=:circle,c=:blue)

n = 4
nplus = 2
Lx = 1.2
Ly = 0.3

κi = [ones(nplus); -ones(n-nplus)]

xd,yd,is_dipole,num_found_in_dipoles = get_dipoles(xi,yi,n,nplus,Lx,Ly)

@test num_found_in_dipoles == 2

## seed_clusters
xi = [-0.549 -0.548 0.160 0.335 -1.772 -1.78]
yi = [0.123 0.16 -0.103 -0.079 0.20 0.220]
xi,yi=xi[:],yi[:]
scatter(xi[1:3],yi[1:3],m=:circle,c=:red)
scatter!(xi[4:6],yi[4:6],m=:circle,c=:blue)
n = 6
nplus = 3
Lx = 1.2
Ly = 0.3
xs,ys,is_seed,seeds_index,num_found_in_seeds = seed_clusters(xi,yi,n,nplus,Lx,Ly)

@test num_found_in_seeds == 4

## grow_plus_clusters

## RCA
# some vortices
xp,yp = xi[1:3],yi[1:3]
xn,yn = xi[4:6],yi[4:6]

# RCA body code
nplus = length(xp)
nminus = length(xn)
n = nplus + nminus
kappai = [ones(nplus); -ones(nminus)]
vortex_index = 1:n



#
#   @info "All vortices are in dipoles: nothing left to do..."
#   PCLUSTERS = struct()
#   NCLUSTERS = struct()
#   DIPOLES.x = x_dipoles
#   DIPOLES.y = y_dipoles
#   DIPOLES.n = size(x_dipoles,1)
#
#   if PLOT
#  PlotDipoles(x_dipoles,y_dipoles,Lx,Ly)
#   end
#   return

# elseif num_found_in_dipoles == n-2 && sum(kappanew == 1) == 1 && sum(kappanew == -1) ==1
#   # NOTE!! This part is a HACK, to fix a BUG currently lurking in MATLAB
#   # 2012b (16/07/2014). RCA fails when there are only two vortices left, each
#   # of different signs (which necessarily must form a dipole). For some reason
#   # MATLAB thinks that the arrays are not the same size when using ==, even
#   # when they are. The problem only occurs when the routine is a function,
#   # so keyboard doesn"t help to debug either. It works fine in the command
#   # line...
#    x_dipoles = [x_dipoles xnew(1) xnew(2)]
#    y_dipoles = [y_dipoles ynew(1) ynew(2)]
#    PCLUSTERS = struct()
#    NCLUSTERS = struct()
#    DIPOLES.x = x_dipoles
#    DIPOLES.y = y_dipoles
#    DIPOLES.n = size(x_dipoles,1)
#
#   # if PLOT
#   #   PlotDipoles(x_dipoles,y_dipoles,Lx,Ly)
#   # end
#
#   @info "All vortices are in dipoles: nothing left to do..."
#   return
# end


## TODO Recursively look for dipoles, until no more can be found
# initial dipole search
x_dipoles,y_dipoles,is_dipole,num_found_in_dipoles = get_dipoles(xi,yi,n,nplus,Lx,Ly)

xnew = @. xi[~Bool(is_dipole)]
ynew = @. yi[~Bool(is_dipole)]
kappanew = @. kappai[~Bool(is_dipole)]
remaining_vortices_index =  @. vortex_index[~Bool(is_dipole)]

np = sum(@. kappanew==1.0)
ntot = length(kappanew)

if num_found_in_dipoles == n
  @info "All vortices are in dipoles. Finished."
end

# while (num_found_in_dipoles >0) && ~isempty(xnew)
  xd,yd,is_dipole2,num_found_in_dipoles = get_dipoles(xnew,ynew,ntot,np,Lx,Ly)

  add_to_dipoles = @. remaining_vortices_index[Bool(is_dipole2)]

  is_dipole = Bool.(is_dipole)
  @. is_dipole[add_to_dipoles] = true

  xnew = @. xi[~is_dipole]
  ynew = @. yi[~is_dipole]
  kappanew = @. kappai[~is_dipole]
  remaining_vortices_index = @. vortex_index[~is_dipole]
  np = sum(@. kappanew==1.0)
  ntot = length(kappanew)
  x_dipoles = [x_dipoles; xd] # ok
  y_dipoles = [y_dipoles; yd] # ok. A little sloppy here ...
# end
