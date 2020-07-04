function PlotDipoles(x_dipoles,y_dipoles,Lx,Ly)
  
    hold on
    
    %Plot Dipoles
    xd = x_dipoles(:,1);
    xd2 = x_dipoles(:,2);
    yd = y_dipoles(:,1);
    yd2 = y_dipoles(:,2);
    plot([xd xd2],[yd yd2],'.g','Markersize',10)
    
    for jj = 1:length(x_dipoles(:,1))
      if (xd2(jj) - xd(jj) > Lx/2) && abs(yd2(jj) - yd(jj)) <= Ly/2
        plot([xd2(jj)-Lx xd(jj)],[ yd2(jj) yd(jj)],'-g','Linewidth',2)
        plot([xd2(jj)  xd(jj)+Lx],[ yd2(jj) yd(jj)],'-g','Linewidth',2)
      elseif (xd2(jj) - xd(jj) < -Lx/2) && abs(yd2(jj) - yd(jj)) <= Ly/2
        plot([xd(jj)-Lx,  xd2(jj)],[ yd(jj), yd2(jj)],'-g','Linewidth',2)
        plot([xd(jj),  xd2(jj)+Lx],[ yd(jj) yd2(jj)],'-g','Linewidth',2)
      elseif (yd2(jj) - yd(jj) > Ly/2) && abs(xd2(jj) - xd(jj)) <= Lx/2
        plot([xd2(jj) xd(jj)],[yd2(jj)-Ly yd(jj)],'-g','Linewidth',2)
        plot([xd2(jj) xd(jj)],[yd2(jj) yd(jj)+Ly],'-g','Linewidth',2)
      elseif (yd2(jj) - yd(jj) < -Ly/2) && abs(xd2(jj) - xd(jj)) <= Lx/2
        plot([xd(jj),  xd2(jj)],[ yd(jj)-Ly, yd2(jj)],'-g','Linewidth',2)
        plot([xd(jj),  xd2(jj)],[ yd(jj) yd2(jj)+Ly],'-g','Linewidth',2)
      elseif (xd2(jj) - xd(jj) > Lx/2) && (yd2(jj) - yd(jj) > Ly/2)
        plot([xd2(jj) xd(jj)+Lx],[yd2(jj) yd(jj)+Ly],'-g','Linewidth',2)
        plot([xd2(jj)-Lx xd(jj)],[yd2(jj)-Ly yd(jj)],'-g','Linewidth',2)
      elseif (xd2(jj) - xd(jj) < -Lx/2) && (yd2(jj) - yd(jj) < -Ly/2)
        plot([xd2(jj) xd(jj)-Lx],[yd2(jj) yd(jj)-Ly],'-g','Linewidth',2)
        plot([xd2(jj)+Lx xd(jj)],[yd2(jj)+Ly yd(jj)],'-g','Linewidth',2)
      elseif (xd2(jj) - xd(jj) < -Lx/2) && (yd2(jj) -yd(jj) > Ly/2)
        plot([xd2(jj) xd(jj)-Lx],[yd2(jj) yd(jj)+Ly],'-g','Linewidth',2)
        plot([xd2(jj)+Lx xd(jj)],[yd2(jj)-Ly yd(jj)],'-g','Linewidth',2)
      elseif (xd2(jj) - xd(jj) > Lx/2) && (yd2(jj) -yd(jj) < -Ly/2)
        plot([xd2(jj) xd(jj)+Lx],[yd2(jj) yd(jj)-Ly],'-g','Linewidth',2)
        plot([xd2(jj)-Lx xd(jj)],[yd2(jj)+Ly yd(jj)],'-g','Linewidth',2)
      else
        plot([ xd(jj) xd2(jj)],[ yd(jj) yd2(jj)],'-g','Linewidth',2)
        
      end
    end
    
    xlim([-Lx/2 Lx/2])
    ylim([-Ly/2 Ly/2])
    box on
    set(gca,'plotboxaspectRatio',[Lx/Ly 1 1])

return