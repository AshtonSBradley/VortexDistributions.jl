using Pkg;Pkg.activate(".")
using Test, Plots, BenchmarkTools, JLD2
using VortexDistributions
gr(transpose=true)

## load test data
loadpath="/Users/abradley/Dropbox/Julia/Vortices - simple results"
f2 = joinpath(loadpath,"nv20-t_200.jld2")

## load good data
@load f2 𝛹200 x y
psi200 = Torus(𝛹200,x,y)
vfound = findvortices(psi200,periodic=true)
p1=heatmap(x,y,angle.(psi200.ψ))
vdata = vortex_array(vfound)
p1

## make a plot function
function plot_vortices!(p,vdata)
    for j in 1:size(vdata)[1]
        vx,vy = vdata[j,1],vdata[j,2]
        scatter!(p1,[vx],[vy],marker=vortex_marker(vdata[j,3]),label=false)
    end
return p1
end

## plot Vortices
plot_vortices!(p1,vdata)

## 
