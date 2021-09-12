using VortexDistributions
using JLD2
using GLMakie
using Interpolations
@load "data/box_rings.jld2" psi_ring1 psi_ring2 psi_ring3 X # make sure julia repl is in base directory 

## Plot iso surface for select ring 
include("util.jl")

psi = psi_ring2;
# plot_iso(psi)

density = abs2.(psi);
pmax = maximum(density);
density = density/pmax;

slice = psi[1:20,1:20,36];
density = abs2.(slice);
pmax = maximum(density);
density = density/pmax
itp = interpolate(density, BSpline(Cubic(Line(OnGrid()))))
heatmap(abs2.(itp), interpolate=true)

scatter!([4.53], [5.4], color="red")
# scatter!([4.6], [5.2], color="red")


itp(4.53, 5.4)

function inBall(v1, v2, R)
    dist = sqrt((v2[1]-v1[1])^2+(v2[2]-v1[2])^2+(v2[3]-v1[3])^2)
    return dist < R
end

## TODO: Check boundary of ball 

