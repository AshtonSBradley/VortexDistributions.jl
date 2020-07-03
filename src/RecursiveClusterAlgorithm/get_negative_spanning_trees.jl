function  NCLUSTERS = GetNegativeSpanningTrees(NCLUSTERS,Lx,Ly)
clusternames = fieldnames(NCLUSTERS);
numclusters = length(clusternames);
for zz = 1:numclusters
  
  positions = eval(['NCLUSTERS.cluster' num2str(zz) '.positions']);
  xi = positions(:,1);
  yi = positions(:,2);
  nc = length(xi);
  %Get the spanning tree
  xdiff= repmat(xi,1,nc) - repmat(xi',nc,1);
  ydiff= repmat(yi,1,nc) - repmat(yi',nc,1);
  rij = sqrt(xdiff.^2 + ydiff.^2);
  rij = tril(rij);
  rij = sparse(rij);
  
  [T,~] = graphminspantree(rij);
  
  %Try wrapping all positions, save the spanning tree which is
  %smallest
  T_old = T;
  [~,I] = sort(xi);
  xi_new = xi;
  yi_new = yi;
  
  %Wrapping in x...
  if ( any(xi<-Lx/4) && any(xi>Lx/4) )
    for ii =  1:length(xi)
      test_xi = xi;
      test_xi(I(1:ii)) = test_xi(I(1:ii)) + Lx;
      xdiff= repmat(test_xi,1,nc) - repmat(test_xi',nc,1);
      ydiff= repmat(yi,1,nc) - repmat(yi',nc,1);
      rij = sqrt(xdiff.^2 + ydiff.^2);
      rij = tril(rij);
      rij = sparse(rij);
      [T,~] = graphminspantree(rij);
      if sum(T(:)) < sum(T_old(:))
        disp('Found smaller spanning tree by wrapping in x')
        T_old = T;
        xi_new = test_xi;
        
      end
    end
  end
  
  %Wrapping in y...
  [~,I] = sort(yi);
  if (any(yi<-Ly/4) && any(yi>Ly/4)) %Changed to any() from sum()>0 (faster) .. 
    for ii =  1:length(yi)
      test_yi = yi;
      test_yi(I(1:ii)) = test_yi(I(1:ii)) + Ly;
      ydiff= repmat(test_yi,1,nc) - repmat(test_yi',nc,1);
      xdiff= repmat(xi_new,1,nc) - repmat(xi_new',nc,1);
      rij = sqrt(xdiff.^2 + ydiff.^2);
      rij = tril(rij);
      rij = sparse(rij);
      [T,~] = graphminspantree(rij);
      if sum(T(:)) < sum(T_old(:))
        disp('Found smaller spanning tree by wrapping in y')
        T_old = T;
        yi_new = test_yi;
      end
    end
  end
  
  %Get centre-of-mass position and radius
  XC = mean(xi_new);
  YC = mean(yi_new);
  Xrad = mean(abs(xi_new-XC));
  Yrad = mean(abs(yi_new-YC));
  Radius = sqrt(Xrad.^2 + Yrad.^2); %#ok
  XC = XC -Lx*(XC > Lx/2);
  YC = YC -Ly*(YC > Ly/2);
  XC = XC +Lx*(XC < -Lx/2);
  YC = YC +Ly*(YC < -Ly/2);
  RC = [XC YC]; %#ok
  
  clear negindex
  [I,J] = find(T_old); %#ok
  eval(['NCLUSTERS.cluster' num2str(zz) '.spanningtreepointers = [I J];']);
  eval(['NCLUSTERS.cluster' num2str(zz) '.RC = RC;']);
  eval(['NCLUSTERS.cluster' num2str(zz) '.Radius = Radius;']);
   
 end
