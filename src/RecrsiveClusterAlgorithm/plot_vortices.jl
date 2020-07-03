function PlotVortices(xi,yi,nplus,nminus,Lx,Ly,R)

n = nplus+nminus;

xp = xi(1:nplus);
yp = yi(1:nplus);
xn = xi(nplus+1:n);
yn = yi(nplus+1:n);


plot(xp,yp,'.r','Markersize',10)
hold on
plot(xn,yn,'.b','Markersize',10)
xlim([-Lx/2 Lx/2])
ylim([-Ly/2 Ly/2])
set(gca,'plotboxaspectratio',[Lx/Ly 1 1])
if nargin == 7
    theta= 0:pi/100:2*pi;
    u = R*cos(theta);
    v = R*sin(theta);
    plot(u,v,'-k','Linewidth',1)
end


