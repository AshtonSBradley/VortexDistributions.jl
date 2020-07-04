
## initialize module
using Pkg;Pkg.activate(".")
using VortexDistributions

## other stuff
using Test, Plots

## set up mock data
N = 10
xi = rand(N)
yi = rand(N)

## make cluster using vortices and tree
vort = PointVortex.(xi,yi,one.(xi))

# default constructor
c0 = Cluster(vort,spanning_tree(xi,yi))

## callable type
c1 = Cluster(vort)
@test c1.vortices |> length == N
# @test c0 == c1

## plot cluster
p1 = plot()
plot_cluster!(p1,c1) # check that the plot method still works

## do some random tests of robustness
N = 30
xi,yi = randn(N),randn(N)
tree = spanning_tree(xi,yi)
c1 = Cluster(PointVortex.(xi,yi,one.(xi)),tree)
p1 = plot(grid=false,axis=false)
plot_cluster!(p1,c1)

xi,yi = randn(N),randn(N)
xi .+= 5
tree = spanning_tree(xi,yi)
c1 = Cluster(PointVortex.(xi,yi,-one.(xi)),tree)
plot_cluster!(p1,c1)

