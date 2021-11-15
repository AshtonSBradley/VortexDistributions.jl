using VortexDistributions, Distributions, Statistics, BenchmarkTools

include("../util.jl")

function euclid(x, y)
    return sqrt((x[1]-y[1])^2 + (x[2]-y[2])^2)
end

##
Nx = 4000; Ny = 4000
Lx = 200; Ly = Lx
x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, Ny+1)[1:end-1];
psi0 = one.(x*y') |> complex
psi = Torus(psi0,x,y)

## 

rand_vorts = rand_pointvortex(2, psi)
# p1 = PointVortex(0., 0., 1)
# p2 = PointVortex(0.001, 0.002, 1)

vortex!(psi, rand_vorts)
# vortex!(psi, [p1, p2])


@time vfound1 = findvortices(psi)
@time vfound2 = naivePlaquette(psi.ψ, [psi.x, psi.y], 0)
@time vfound3 = remove_vortices_edge(findvortices_jumps(psi), psi)
@time vfound4 = remove_vortices_edge(findvortices_grid(psi), psi)

# sort!(rand_vorts, by = x -> x.xv)
sort!(vfound1, by = x -> x.xv)
sort!(vfound2, by = x -> x.xv)
sort!(vfound3, by = x -> x.xv)
sort!(vfound4, by = x -> x.xv)

vf1diff = [euclid([rand_vorts[i].xv, rand_vorts[i].yv], [vfound1[i].xv, vfound1[i].yv]) for i in 1:length(rand_vorts)]
vf2diff = [euclid([rand_vorts[i].xv, rand_vorts[i].yv], [vfound2[i].xv, vfound2[i].yv]) for i in 1:length(rand_vorts)]
vf2_3diff = [euclid([vfound2[i].xv, vfound2[i].yv], [vfound3[i].xv, vfound3[i].yv]) for i in 1:length(rand_vorts)]


vf1avg = mean(vf1diff)
vf2avg = mean(vf2diff)
vf2_3avg = mean(vf2_3diff)

naive = benchmark2Dnaive(2)

naive = @benchmark naivePlaquette(psi.ψ, [psi.x, psi.y], 0)
optimised = @benchmark remove_vortices_edge(findvortices_grid(psi), psi)

benchmarks = benchmark2Dnaive2(2)

using JLD2

@save "data/naiveBenchmarks.jld2" benchmarks;