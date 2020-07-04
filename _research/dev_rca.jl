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
@test c0.vortices[1] == c1.vortices[1]
@test c0.tree[1] == c1.tree[1]

## plot cluster
p1 = plot()
plot_cluster!(p1,c1) # check that the plot method still works

## do some random tests of robustness
N = 20
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

## vector of clusters of different sizes
xi,yi = randn(2N),randn(2N)
xi .+= 10
tree = spanning_tree(xi,yi)
c2 = Cluster(PointVortex.(xi,yi,one.(xi)),tree)
cvec = [c0;c1;c2]

xi,yi = randn(3N),randn(3N)
xi .+= 20
tree = spanning_tree(xi,yi)
c3 = Cluster(PointVortex.(xi,yi,-one.(xi)),tree)
push!(cvec,c3)

## cvec
Base.size(c::Cluster) = length(c.vortices)
cvec[3] |> size
size(cvec)

## grow_plus_clusters
