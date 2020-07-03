## initialize module
using Pkg;Pkg.activate(".")
# include("../src/RecursiveClusterAlgorithm.jl")
# using .RecursiveClusterAlgorithm

using VortexDistributions
using Test, Plots

## set up mock data
N = 10
xi = rand(N)
yi = rand(N)
tree = spanning_tree(xi,yi)

## make cluster using vortices and tree
vort = PointVortex.(xi,yi,one.(xi))
c1 = Cluster(vort,tree)
@test c1.vortices |> length == N

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

