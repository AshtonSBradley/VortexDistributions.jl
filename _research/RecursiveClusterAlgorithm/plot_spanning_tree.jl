function PlotSpanningTree(xi,yi,clustersign,Lx,Ly,treepointers,fighandle) 

I = treepointers(:,1);
J = treepointers(:,2);

x1 = xi(I);
x2 = xi(J);
y1 = yi(I);
y2 = yi(J);

if clustersign == 1
    c = [1 0 0];
elseif clustersign == -1
    c = [0 0 1];
end

lw = 1.5;
%figure(fighandle)
hold on

for ii = 1:length(I)
    if (x2(ii) - x1(ii) < -Lx/2) && ( abs(y2(ii) - y1(ii)) < Ly/2)
        plot([x1(ii)-Lx,  x2(ii)],[ y1(ii), y2(ii)],'-','Linewidth',lw,'Color',c)
        plot([x1(ii),  x2(ii)+Lx],[ y1(ii) y2(ii)],'-','Linewidth',lw,'Color',c)
    elseif (x2(ii) - x1(ii) > Lx/2) && (abs(y2(ii) - y1(ii)) < Ly/2)
        plot([x2(ii)-Lx,  x1(ii)],[ y2(ii), y1(ii)],'-','Linewidth',lw,'Color',c)
        plot([x2(ii),  x1(ii)+Lx],[ y2(ii) y1(ii)],'-','Linewidth',lw,'Color',c)
    elseif (y2(ii) - y1(ii) < -Ly/2) && (abs(x2(ii) - x1(ii)) < Lx/2)
        plot([x1(ii) x2(ii)],[y1(ii) y2(ii)+Ly],'-','Linewidth',lw,'Color',c)
        plot([x2(ii) x1(ii)],[y2(ii) y1(ii)-Ly],'-','Linewidth',lw,'Color',c)
    elseif (y2(ii) - y1(ii) > Ly/2) && (abs(x2(ii) - x1(ii)) < Lx/2)
        plot([x1(ii) x2(ii)],[y1(ii) y2(ii)-Ly],'-','Linewidth',lw,'Color',c)
        plot([x2(ii) x1(ii)],[y2(ii) y1(ii)+Ly],'-','Linewidth',lw,'Color',c)
    elseif (y2(ii) - y1(ii) < -Ly/2) && (x2(ii) - x1(ii) < -Lx/2)
        plot([x1(ii)-Lx,  x2(ii)],[ y1(ii)-Ly, y2(ii)],'-','Linewidth',lw,'Color',c)
        plot([x1(ii),  x2(ii)+Lx],[ y1(ii) y2(ii)+Ly],'-','Linewidth',lw,'Color',c)
    elseif (x2(ii) - x1(ii) > Lx/2) && (y2(ii) - y1(ii) < -Ly/2)
        plot([x2(ii)-Lx,  x1(ii)],[ y2(ii)+Ly, y1(ii)],'-','Linewidth',lw,'Color',c)
        plot([x2(ii),  x1(ii)+Lx],[ y2(ii) y1(ii)-Ly],'-','Linewidth',lw,'Color',c)
    elseif (x2(ii) - x1(ii) < -Lx/2) && (y2(ii) - y1(ii) > Ly/2)
        plot([x2(ii)+Lx,  x1(ii)],[ y2(ii)-Ly, y1(ii)],'-','Linewidth',lw,'Color',c)
        plot([x2(ii),  x1(ii)-Lx],[ y2(ii) y1(ii)+Ly],'-','Linewidth',lw,'Color',c)
    elseif (x2(ii) - x1(ii) < -Lx/2) && (y2(ii)-y1(ii) <-Ly/2)
        plot([x2(ii)+Lx,  x1(ii)],[ y2(ii)+Ly, y1(ii)],'-','Linewidth',lw,'Color',c)
        plot([x2(ii),  x1(ii)-Lx],[ y2(ii) y1(ii)-Ly],'-','Linewidth',lw,'Color',c)
    elseif (x2(ii) - x1(ii) > Lx/2) && (y2(ii)-y1(ii) > Ly/2)
        plot([x2(ii),  x1(ii)+Lx],[ y2(ii), y1(ii)+Ly],'-','Linewidth',lw,'Color',c)
        plot([x2(ii)-Lx,  x1(ii)],[ y2(ii)-Ly y1(ii)],'-','Linewidth',lw,'Color',c)
      
    else
        plot([x1(ii) x2(ii)],[y1(ii) y2(ii)],'-','Linewidth',lw,'Color',c)
    end
end

return