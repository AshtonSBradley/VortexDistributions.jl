%function FindSmallestSpanningTree(numclusters,PCLUSTERS) %#ok
numclusters= cluster_id_number;
for zz = 1:numclusters
    positions = eval(['PCLUSTERS.cluster' num2str(zz) '.positions']);
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
    [I,J] = find(T);
    
    figure(1)
    hold on
    for ii = 1:length(I)
        plot([xi(I(ii)) xi(J(ii))],[yi(I(ii)) yi(J(ii))],'-r','Linewidth',2)
    end
    
    %Try wrapping all negative positions, save the spanning tree which is
    %smallest
    T_old = T;
    
    %Wrapping in x...
    if ( (sum(xi<-Lx/4)>0) && (sum(xi>Lx/4)>0) )
      negindex = (xi<0);
      posindex = (xi>0);
      numnegs = sum(negindex);
      numpos = sum(posindex);
      posindex = find(posindex);
      negindex = find(negindex);
      keyboard
      if (numnegs <= numpos)
        for ii =  1:numnegs
          test_xi = xi;
          test_xi(negindex(1:ii)) = test_xi(negindex(1:ii)) + Lx;
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
            yi_new = yi;
          end
        end
      elseif (numnegs > numpos)
        
        for ii =  1:numpos
          test_xi = xi;
          test_xi(posindex(1:ii)) = test_xi(posindex(1:ii)) - Lx;
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
            yi_new = yi;
          end
        end
      end
    else
      xi_new = xi;
      yi_new = yi;
    end
    
    %Wrapping in y    
    if ((sum(yi<-Ly/4)>0) && (sum(yi>Ly/4)>0))
      keyboard
      negindex = (yi<0);
      posindex = (yi>0);
      numnegs = sum(negindex);
      numpos = sum(posindex);
      negindex = find(negindex);
      posindex = find(posindex);
      if (numnegs <= numpos)
        for ii =  1:numnegs
          test_yi = yi;
          test_yi(negindex(1:ii)) = test_yi(negindex(1:ii)) + Ly;
          ydiff= repmat(test_yi,1,nc) - repmat(test_yi',nc,1);
          xdiff= repmat(xi,1,nc) - repmat(xi',nc,1);
          rij = sqrt(xdiff.^2 + ydiff.^2);
          rij = tril(rij);
          rij = sparse(rij);
          [T,~] = graphminspantree(rij);
          if sum(T(:)) < sum(T_old(:))
            disp('Found smaller spanning tree by wrapping in y')
            T_old = T;
            yi_new = test_yi;
            xi_new = xi;
          end
        end
      elseif (numnegs > numpos)
        for ii =  1:numpos
          test_yi = yi;
          test_yi(posindex(1:ii)) = test_yi(posindex(1:ii)) - Ly;
          ydiff= repmat(test_yi,1,nc) - repmat(test_yi',nc,1);
          xdiff= repmat(xi,1,nc) - repmat(xi',nc,1);
          rij = sqrt(xdiff.^2 + ydiff.^2);
          rij = tril(rij);
          rij = sparse(rij);
          [T,~] = graphminspantree(rij);
          if sum(T(:)) < sum(T_old(:))
            disp('Found smaller spanning tree by wrapping in y')
            T_old = T;
            yi_new = test_yi;
            xi_new = xi;
          end
        end
      end
    else
      xi_new = xi;
      yi_new = yi;
    end
    
    clear negindex 
    [I,J] = find(T_old);
   eval(['PCLUSTERS.cluster' num2str(zz) '.spanningtreepointers = [I J];']);
  
   %Rewrap
    xi_new(xi_new > Lx/2) = xi_new(xi_new > Lx/2) - Lx;
    yi_new(yi_new > Ly/2) = yi_new(yi_new > Ly/2) - Ly;
    xi_new(xi_new < Lx/2) = xi_new(xi_new < Lx/2) + Lx;
    yi_new(yi_new < Ly/2) = yi_new(yi_new < Ly/2) + Ly;
    
        x1 = xi(I);
    x2 = xi(J);
    y1 = yi(I);
    y2 = yi(J);
         
    
    figure(2)
    hold on
    for ii = 1:length(I)
        if (x2(ii) - x1(ii) < -Lx/2) && ( abs(y2(ii) - y1(ii)) < Ly/2)
            plot([x1(ii)-Lx,  x2(ii)],[ y1(ii), y2(ii)],'-r','Linewidth',2)
            plot([x1(ii),  x2(ii)+Lx],[ y1(ii) y2(ii)],'-r','Linewidth',2)
        elseif (x2(ii) - x1(ii) > Lx/2) && (abs(y2(ii) - y1(ii)) < Ly/2)
            plot([x2(ii)-Lx,  x1(ii)],[ y2(ii), y1(ii)],'-r','Linewidth',2)
            plot([x2(ii),  x1(ii)+Lx],[ y2(ii) y1(ii)],'-r','Linewidth',2)
        elseif (y2(ii) - y1(ii) < -Ly/2) && (abs(x2(ii) - x1(ii)) < Lx/2)
            plot([x1(ii) x2(ii)],[y1(ii) y2(ii)+Ly],'-r','Linewidth',2)
            plot([x2(ii) x1(ii)],[y2(ii) y1(ii)-Ly],'-r','Linewidth',2)
        elseif (y2(ii) - y1(ii) > Ly/2) && (abs(x2(ii) - x1(ii)) < Lx/2)
            plot([x1(ii) x2(ii)],[y1(ii) y2(ii)-Ly],'-r','Linewidth',2)
            plot([x2(ii) x1(ii)],[y2(ii) y1(ii)+Ly],'-r','Linewidth',2)
        elseif (y2(ii) - y1(ii) < -Ly/2) && (x2(ii) - x1(ii) < -Lx/2)
            plot([x1(ii)-Lx,  x2(ii)],[ y1(ii)-Ly, y2(ii)],'-r','Linewidth',2)
            plot([x1(ii),  x2(ii)+Lx],[ y1(ii) y2(ii)+Ly],'-r','Linewidth',2)
        else
            plot([x1(ii) x2(ii)],[y1(ii) y2(ii)],'-r','Linewidth',2)
        end
    end
    
end
clear xi_new yi_new x1 x2
xlim([-Lx/2 Lx/2])
ylim([-Ly/2 Ly/2])
axis('square')













