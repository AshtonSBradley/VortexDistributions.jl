function grow_plus_clusters(xnew,ynew,kappanew,seeds_index,seednum,Lx,Ly)

xseeds = xnew[seeds_index]
yseeds = ynew[seeds_index]
size_check = size(xseeds)
#!! For some reason, on the last iteration the dimensionality changes!? !!#
if  size_check[2] == 1
    xseeds = xseeds'
    yseeds = yseeds'
end

xc = xseeds[seednum,:]
yc = yseeds[seednum,:]
nc = 2 # initial cluster charge is 2

# define current vortex indices associated with the cluster
cluster_index = seeds_index[seednum,:]
num_added_to_cluster = nc # seeds are added to the cluster

xp = @. xnew[kappanew == 1.0]
yp = @. ynew[kappanew == 1.0]
xn = @. xnew[kappanew == -1.0]
yn = @. ynew[kappanew == -1.0]
failed_attempts = 0

#TODO mapped to at least here

while (num_added_to_cluster > 0 && failed_attempts < 10)
    # ignore vortices already in the cluster
        xtemp = xp
        ytemp = yp
        for j in 1:nc
            Ij = ( (xc[j] == xp) && (yc[j] == yp) )
            xtemp[Ij] = NaN
            ytemp[Ij] = NaN
            # clear I

            # Work out the distances to all other vortices
            xdiff = @. xtemp - xc'
            ydiff = @. ytemp - yc'
        end

    # incorporate periodicity
    xdiff,ydiff = periodic_distances(xdiff,ydiff,Lx,Ly)
    rij = hypot.(xdiff,ydiff)

    # Nearest same-sign vortex to each vortex, numbered by x
    rmin,candidate_index = minimum(rij)
    # candidate_index = unique(candidate_index)

   ## Find the nearest opposites to the cluster
    xdiff3 = @. xn - xc'
    ydiff3 = @. yn - yc'
    rij3 = hypot.(xdiff3,ydiff3)

    r_min3,_ = minimum(rij3)

   # candidate_index = candidate_index(rmin < r_min3)

    ## Find the nearest opposite signed vortex of each candidate
    number_of_candidates = length(candidate_index)

    x_candidates = xp[candidate_index]
    y_candidates = yp[candidate_index]
    xdiff2 = @. xn - x_candidates'
    ydiff2 = @. yn - y_candidates'
    xdiff2,ydiff2 = periodic_distances(xdiff2,ydiff2,Lx,Ly)

    rij2 = hypot.(xdiff2,ydiff2)
    rmin_opposite,_ = minimum(rij2)

    add_to_cluster = ( (rmin < rmin_opposite) && (rmin < r_min3) )
    num_added_to_cluster = length(unique(candidate_index[add_to_cluster]))
    cluster_index = [cluster_index unique(candidate_index[add_to_cluster])] ##ok

    xc = xnew[cluster_index]'
    yc = ynew[cluster_index]'
    nc = length(xc)
#     disp(['found ' num2str(number_of_candidates) ' candidates'])
#     disp(['added ' num2str(num_added_to_cluster) ' candidates to cluster'])
end
xc = xc'
yc = yc'
cluster_index = cluster_index'
return xc,yc,cluster_index
end
