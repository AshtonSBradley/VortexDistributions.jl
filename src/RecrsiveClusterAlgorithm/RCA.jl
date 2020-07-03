## Written by Matt Reeves 2013
## RCA -- Recursive Cluster Algorithm
#    The RCA takes a given distribution of quantum vortices and decomposes it
#    into dipoles,clusters and free vortices, as described in Reeves et al, PRL
#    2012
#
## Inputs:
#   rp - vector of positive vortex positions [xp1 yp1 xp2 yp2 ... xpN ypN]
#   rn - vector of positive vortex positions [x1 y1 x2 y2 ... xnN ynN]
#   Lx/Ly - spatial domain size in healing lengths \xi = \hbar/\sqrt(nmg)
#   PLOT - optional logical that specifies if the RCA spanning tree
#   is to be plotted
#
## Outputs
#   PCLUSTERS, NCLUSTERS - structs containing the charge, radius, centre of
#   mass positions of each cluster...

function RCA(rp,rn,Lx,Ly,PLOT)

## TODO -- FIX BUG WHERE THINGS FAIL IF EVERYTHING IS IN A DIPOLE!!!!!!!!!!
##

## Unpack positions, get vortex numbers, etc.
#rp = dlmread([rootfile "PositivePositions.txt"])
#rn = dlmread([rootfile "NegativePositions.txt"])
xp = rp[:,1]
yp = rp[:,2]
xn = rn[:,1]
yn = rn[:,2]
xi = [xp; xn]
yi = [yp; yn]

nplus = length(xp)
nminus = length(xn)
n = nplus + nminus
kappai = [ones(nplus); -ones(nminus)]
vortex_index = 1:n
## Recursively Remove Dipoles (Mutual, opposite sign nearest neighbours)

#Initial dipole search
x_dipoles,y_dipoles,is_dipole,num_found_in_dipoles = get_dipoles(xi,yi,n,nplus,Lx,Ly)
xnew = @. xi[~is_dipole]
ynew = @. yi[~is_dipole]
kappanew = @. kappai[~is_dipole]
remaining_vortices_index =  vortex_index[~is_dipole]
np = sum(kappanew==1)
ntot = length(kappanew)

if num_found_in_dipoles == n

  disp("All vortices are in dipoles: nothing left to do...")
  PCLUSTERS = struct()
  NCLUSTERS = struct()
  DIPOLES.x = x_dipoles
  DIPOLES.y = y_dipoles
  DIPOLES.n = size(x_dipoles,1)

  if PLOT
 PlotDipoles(x_dipoles,y_dipoles,Lx,Ly)
  end
  return

elseif num_found_in_dipoles == n-2 && sum(kappanew == 1) == 1 && sum(kappanew == -1) ==1
  # NOTE!! This part is a HACK, to fix a BUG currently lurking in MATLAB
  # 2012b (16/07/2014). RCA fails when there are only two vortices left, each
  # of different signs (which necessarily must form a dipole). For some reason
  # MATLAB thinks that the arrays are not the same size when using ==, even
  # when they are. The problem only occurs when the routine is a function,
  # so keyboard doesn"t help to debug either. It works fine in the command
  # line...
   x_dipoles = [x_dipoles xnew(1) xnew(2)]
   y_dipoles = [y_dipoles ynew(1) ynew(2)]
   PCLUSTERS = struct()
   NCLUSTERS = struct()
   DIPOLES.x = x_dipoles
   DIPOLES.y = y_dipoles
   DIPOLES.n = size(x_dipoles,1)

  if PLOT
    PlotDipoles(x_dipoles,y_dipoles,Lx,Ly)
  end

  disp("All vortices are in dipoles: nothing left to do...")
  return
end

#Recursively look for dipoles, until no more can be found
while num_found_in_dipoles >0 && ~isempty(xnew)
  [xd,yd,is_dipole2,num_found_in_dipoles] = GetDipoles(xnew,ynew,ntot,np,Lx,Ly)
  add_to_dipoles = remaining_vortices_index(is_dipole2)
  is_dipole(add_to_dipoles) = true
  xnew = xi(~is_dipole)
  ynew = yi(~is_dipole)
  kappanew = kappai(~is_dipole)
  remaining_vortices_index = vortex_index(~is_dipole)
  np = sum(kappanew==1)
  ntot = length(kappanew)
  x_dipoles = [x_dipoles xd]##ok
  y_dipoles = [y_dipoles yd]##ok  #A little sloppy here...
end

DIPOLES.x = x_dipoles
DIPOLES.y = y_dipoles
DIPOLES.n = size(x_dipoles,1)

if isempty(xnew)
  PCLUSTERS = struct()
  NCLUSTERS = struct()

  if PLOT
   PlotDipoles(x_dipoles,y_dipoles,Lx,Ly)

  end

  return
end

## Find Cluster Seeds (Mutual, same sign nearest-neighbours)
if length(xnew) == 1
  disp("Only one vortex left nothing left to do...")
  PCLUSTERS = struct()
  NCLUSTERS = struct()

  if PLOT
    PlotDipoles(x_dipoles,y_dipoles,Lx,Ly)
  end

  return
end
[xseeds,yseeds,is_seed,seeds_index,num_found_in_seeds] = SeedClusters(xnew,ynew,ntot,np,Lx,Ly)
kappa_seeds = kappanew(is_seed)

if num_found_in_seeds >0
  #Remove doubles in seeds
  #seeds_index_old = seeds_index     #check that the below does not remove any indices incorrectly
  counter = 0
  check = 1
  while sum(check)~=0 && counter < length(xseeds(:,1)) #DOING THIS WITH A FOR LOOP WOULD BE BETTER
    counter = counter+1
    check =  ( (seeds_index(counter,1) == seeds_index(:,2)) & (seeds_index(counter,2) == seeds_index(:,1)) )
    I = find(check==0)
    xseeds = xseeds(I,:)
    yseeds = yseeds(I,:)
    kappa_seeds = kappa_seeds(I)
    seeds_index = seeds_index(I,:)
  end
  clear I
  ## Grow Plus Clusters From seeds
  if sum(kappanew(is_seed) ==1) >0
    disp("Growing Positive Clusters...")
    cluster_id_number = 1
    seednum = 1 #note that seednum is always 1
    posind = kappa_seeds == 1
    pos_seeds_index = [seeds_index(posind,1) seeds_index(posind,2)]
    [xc,yc,cluster_index] = GrowPlusClusters2(xnew,ynew,kappanew,pos_seeds_index,seednum,Lx,Ly)##ok

    eval(["PCLUSTERS.cluster" num2str(cluster_id_number) ".positions = [xc yc]"])
    eval(["PCLUSTERS.cluster" num2str(cluster_id_number) ".indices = cluster_index"])

    # check which seeds have been used
    checkvec = [cluster_index pos_seeds_index(:)]
    [duplicates,duplicate_indices] = hist(checkvec,1:max(checkvec))
    duplicate_indices = duplicate_indices(duplicates == 2)
    if sum(duplicates >2) > 0
      error("The same seed has been added to the cluster more than once?!?!")
    end

    #remove used seeds from consideration
    [~,I] = setdiff(pos_seeds_index(:,1),duplicate_indices)
    seeds_index_new = [pos_seeds_index(I,1) pos_seeds_index(I,2)]

    num_used_seeds = length(pos_seeds_index(:,1)) - length(seeds_index_new(:,1))
    num_remaining_seeds = length(seeds_index_new(:,1))
    disp(["Used " num2str(num_used_seeds) " seeds in last cluster growth: " ...
      num2str(num_remaining_seeds) " seeds remaining"])

    while num_remaining_seeds > 0
      #Grow a new cluster, using the first seed of the new seed vector
      [xc,yc,cluster_index] = ...
        GrowPlusClusters2(xnew,ynew,kappanew,seeds_index_new,seednum,Lx,Ly) ##ok


      # !! NOTE !! - A multiple merge has never been a problem but did
      # recently occur, at least once, for a VERY large number of
      #vortices (8192). I suspect in this case, either:
      #       a)  Merging should be moved to the end of the growing procedure
      #       b) a check should be performed at the end of the growing procedure. MTR 05/03/2015.
      num_merges = 0
      for ii = 1:cluster_id_number
        check1 = eval(["length(unique([PCLUSTERS.cluster" num2str(ii) ".indices  cluster_index]))"])
        check2 = eval(["length([PCLUSTERS.cluster" num2str(ii) ".indices cluster_index])"])
        if check1 ~= check2
          disp(["Cluster " num2str(ii) " and current cluster (" ...
            num2str(cluster_id_number+1) ") are not unique: merging..."])
          num_merges = num_merges+1
          if num_merges > 1
            disp("WARNING: a cluster has been merged with two different clusters -- please read comment at line 184")
          end
          cluster_index = eval(["unique([cluster_index PCLUSTERS.cluster" num2str(ii) ".indices])"])
          xc = xnew(cluster_index)##ok
          yc = ynew(cluster_index)##ok
          eval(["PCLUSTERS.cluster" num2str(ii) ".positions = [xc yc]"])
          eval(["PCLUSTERS.cluster" num2str(ii) ".indices = cluster_index"])
        end
        clear check1 check2
      end

      if num_merges == 0
        cluster_id_number = cluster_id_number +1
        disp(["New cluster formed: " num2str(cluster_id_number)])
        eval(["PCLUSTERS.cluster" num2str(cluster_id_number) ".positions = [xc yc]"])
        eval(["PCLUSTERS.cluster" num2str(cluster_id_number) ".indices = cluster_index"])
      end

      # check which seeds have been used
      checkvec = [cluster_index pos_seeds_index(:)]
      [duplicates,duplicate_indices] = hist(checkvec,1:max(checkvec))
      duplicate_indices = duplicate_indices(duplicates == 2)
      if sum(duplicates >2) > 0
        error("The same seed has been added to the cluster more than once?!?!")
      end

      # remove used seeds from consideration
      [~,I] = setdiff(seeds_index_new(:,1),duplicate_indices)
      seeds_index_new = [seeds_index_new(I,1) seeds_index_new(I,2)]

      num_used_seeds = length(pos_seeds_index(:,1)) - length(seeds_index_new(:,1))
      num_remaining_seeds = length(seeds_index_new(:,1))
      disp(["Used " num2str(num_used_seeds) " seeds: " ...
        num2str(num_remaining_seeds) " seeds remaining"])

    end

    PCLUSTERS = GetPositiveSpanningTrees(PCLUSTERS,Lx,Ly)##ok
  else
    PCLUSTERS = struct()
  end
  ## Grow Minus Clusters From seeds
  if sum(kappanew(is_seed) ==-1) >0
    disp("Growing Negative Clusters...")
    cluster_id_number = 1
    seednum = 1
    negind = kappa_seeds == -1
    neg_seeds_index = [seeds_index(negind,1) seeds_index(negind,2)]
    [xc,yc,cluster_index] = ...
      GrowMinusClusters2(xnew,ynew,kappanew,neg_seeds_index,seednum,Lx,Ly)   ##ok

    eval(["NCLUSTERS.cluster" num2str(cluster_id_number) ".positions = [xc yc]"])
    eval(["NCLUSTERS.cluster" num2str(cluster_id_number) ".indices = cluster_index"])

    # check which seeds have been used
    checkvec = [cluster_index neg_seeds_index(:)]
    [duplicates,duplicate_indices] = hist(checkvec,1:max(checkvec))
    duplicate_indices = duplicate_indices(duplicates == 2)
    if sum(duplicates >2) > 0
      error("The same seed has been added to the cluster more than once?!?!")
    end

    #remove used seeds from consideration
    [~,I] = setdiff(neg_seeds_index(:,1),duplicate_indices)
    seeds_index_new = [neg_seeds_index(I,1) neg_seeds_index(I,2)]

    num_used_seeds = length(neg_seeds_index(:,1)) - length(seeds_index_new(:,1))
    num_remaining_seeds = length(seeds_index_new(:,1))
    disp(["Used " num2str(num_used_seeds) " seeds in last cluster growth: " ...
      num2str(num_remaining_seeds) " seeds remaining"])

    while num_remaining_seeds > 0
      #Grow a new cluster, using the first seed of the new seed vector
      [xc,yc,cluster_index] = ...
        GrowMinusClusters2(xnew,ynew,kappanew,seeds_index_new,seednum,Lx,Ly) ##ok
      num_merges = 0

      for ii = 1:cluster_id_number
        check1 = eval(["length(unique([NCLUSTERS.cluster" num2str(ii) ".indices  cluster_index]))"])
        check2 = eval(["length([NCLUSTERS.cluster" num2str(ii) ".indices cluster_index])"])
        if check1 ~= check2
          disp(["Cluster " num2str(ii) " and current cluster (" num2str(cluster_id_number+1) ") are not unique: merging..."])
          num_merges = num_merges+1
          if num_merges > 1
           disp("WARNING: a cluster has been merged with two different clusters -- please read comment at line 184")
          end
          cluster_index = eval(["unique([cluster_index NCLUSTERS.cluster" num2str(ii) ".indices])"])
          xc = xnew(cluster_index)##ok
          yc = ynew(cluster_index)##ok
          eval(["NCLUSTERS.cluster" num2str(ii) ".positions = [xc yc]"])
          eval(["NCLUSTERS.cluster" num2str(ii) ".indices = cluster_index"])
        end
        clear check1 check2
      end

      if num_merges == 0
        cluster_id_number = cluster_id_number +1
        disp(["New cluster formed: " num2str(cluster_id_number)])
        eval(["NCLUSTERS.cluster" num2str(cluster_id_number) ".positions = [xc yc]"])
        eval(["NCLUSTERS.cluster" num2str(cluster_id_number) ".indices = cluster_index"])
      end

      # check which seeds have been used
      checkvec = [cluster_index neg_seeds_index(:)]
      [duplicates,duplicate_indices] = hist(checkvec,1:max(checkvec))
      duplicate_indices = duplicate_indices(duplicates == 2)
      if sum(duplicates >2) > 0
        error("The same seed has been added to the cluster more than once?!?!")
      end

      ## Remove used seeds from consideration
      [~,I] = setdiff(seeds_index_new(:,1),duplicate_indices)
      seeds_index_new = [seeds_index_new(I,1) seeds_index_new(I,2)]
      num_used_seeds = length(neg_seeds_index(:,1)) - length(seeds_index_new(:,1))
      num_remaining_seeds = length(seeds_index_new(:,1))

      disp(["Used " num2str(num_used_seeds) " seeds: " ...
        num2str(num_remaining_seeds) " seeds remaining"])

    end

    NCLUSTERS = GetNegativeSpanningTrees(NCLUSTERS,Lx,Ly)##ok
  else NCLUSTERS = struct()
  end
else
  disp("No Clusters found, all vortices are in dipoles or free...")
  PCLUSTERS = struct()
  NCLUSTERS = struct()
end

## Plotting to check
if PLOT
  # Plot all the vortices
  #figure(gcf)
  PlotVortices(xi,yi,nplus,nminus,512,512)
  hold on

  #Plus Clusters
  clusternames = fieldnames(PCLUSTERS)
  numclusters = length(clusternames)
  for zz = 1:numclusters
    mystring = ["PCLUSTERS.cluster" num2str(zz)]
    positions = eval([mystring ".positions"])
    pointers = eval([mystring  ".spanningtreepointers"])
    PlotSpanningTree(positions(:,1),positions(:,2),1,Lx,Ly,pointers,gcf)
  end
  xlim([-Lx/2 Lx/2])
  ylim([-Lx/2 Lx/2])

  #Minus Clusters
  clusternames = fieldnames(NCLUSTERS)
  numclusters = length(clusternames)
  for zz = 1:numclusters
    mystring = ["NCLUSTERS.cluster" num2str(zz)]
    positions = eval([mystring ".positions"])
    pointers = eval([mystring  ".spanningtreepointers"])
    PlotSpanningTree(positions(:,1),positions(:,2),-1,Lx,Ly,pointers,gcf)
  end

  PlotDipoles(x_dipoles,y_dipoles,Lx,Ly)
end

return PCLUSTERS, NCLUSTERS, DIPOLES
end
