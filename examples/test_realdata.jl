using Pkg
pkg"activate ."

using Plots, JLD2, Test, Revise, VortexDistributions

@load "./examples/one_frame.jld2"

heatmap(x,y,abs2.(ψ1'))

nt,np,nn,vortices = findvortices(ψ1,x,y)
