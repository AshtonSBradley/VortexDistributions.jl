using Plots, VortexDistributions

@load "./examples/one_frame.jld2"

heatmap(x,y,abs2.(ψ1'))

psi = Torus(ψ1,x,y)
vortices = findvortices(psi)

# zoom window
x0,y0 = 70,200
Lx,Ly = 100,100

# simple cartesian mask
xm = x[@. abs(x-x0)<Lx/2]
ym = y[@. abs(y-y0)<Ly/2]
psim = @. ψ1[abs(x-x0)<Lx/2,abs(y-y0)<Ly/2]

heatmap(xm,ym,angle.(psim),transpose=true)
psi2 = Torus(psim,xm,ym)
vort = findvortices(psi2)

scatter!([vort[1].xv],[vort[1].yv],alpha=0.4)
scatter!([vort[2].xv],[vort[2].yv],alpha=0.4)
