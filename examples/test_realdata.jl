using Pkg
pkg"activate ."

using VortexDistributions, Plots, JLD2, Revise

@load "one_frame.jld2"

nt,np,nn,vortices = findvortices(ψ1,x,y)
