using VortexDistributions, Distributions, Statistics, BenchmarkTools

include("../util.jl")

##
N = 2^4
Lx = 200; Ly = Lx
x = LinRange(-Lx/2,Lx/2, N+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, N+1)[1:end-1];
psi0 = one.(x*y') |> complex
psi4 = Torus(psi0,x,y)

##
N=2^5
x = LinRange(-Lx/2,Lx/2, N+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, N+1)[1:end-1];
psi0 = one.(x*y') |> complex
psi5 = Torus(psi0,x,y)

##
N=2^6
x = LinRange(-Lx/2,Lx/2, N+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, N+1)[1:end-1];
psi0 = one.(x*y') |> complex
psi6 = Torus(psi0,x,y)

##
N=2^7
x = LinRange(-Lx/2,Lx/2, N+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, N+1)[1:end-1];
psi0 = one.(x*y') |> complex
psi7 = Torus(psi0,x,y)

##
N=2^8
x = LinRange(-Lx/2,Lx/2, N+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, N+1)[1:end-1];
psi0 = one.(x*y') |> complex
psi8 = Torus(psi0,x,y)

##
N=2^9
x = LinRange(-Lx/2,Lx/2, N+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, N+1)[1:end-1];
psi0 = one.(x*y') |> complex
psi9 = Torus(psi0,x,y)

##
N=2^10
x = LinRange(-Lx/2,Lx/2, N+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, N+1)[1:end-1];
psi0 = one.(x*y') |> complex
psi10 = Torus(psi0,x,y)

##
N=2^11
x = LinRange(-Lx/2,Lx/2, N+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, N+1)[1:end-1];
psi0 = one.(x*y') |> complex
psi11 = Torus(psi0,x,y)

##
N=2^12
x = LinRange(-Lx/2,Lx/2, N+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, N+1)[1:end-1];
psi0 = one.(x*y') |> complex
psi12 = Torus(psi0,x,y)


##
rand_vorts = rand_pointvortex(2, psi)

vortex!(psi4, rand_vorts)
vortex!(psi5, rand_vorts)
vortex!(psi6, rand_vorts)
vortex!(psi7, rand_vorts)
vortex!(psi8, rand_vorts)
vortex!(psi9, rand_vorts)
vortex!(psi10, rand_vorts)
vortex!(psi11, rand_vorts)
vortex!(psi12, rand_vorts)



bench4 = @benchmark naivePlaquette(psi4.ψ, [psi4.x, psi4.y], 0)
bench5 = @benchmark naivePlaquette(psi5.ψ, [psi5.x, psi5.y], 0)
bench6 = @benchmark naivePlaquette(psi6.ψ, [psi6.x, psi6.y], 0)
bench7 = @benchmark naivePlaquette(psi7.ψ, [psi7.x, psi7.y], 0)
bench8 = @benchmark naivePlaquette(psi8.ψ, [psi8.x, psi8.y], 0)
bench9 = @benchmark naivePlaquette(psi9.ψ, [psi9.x, psi9.y], 0)
bench10 = @benchmark naivePlaquette(psi10.ψ, [psi10.x, psi10.y], 0)
bench11 = @benchmark naivePlaquette(psi11.ψ, [psi11.x, psi11.y], 0)
bench12 = @benchmark naivePlaquette(psi12.ψ, [psi12.x, psi12.y], 0)

bench_all = [bench4, bench5, bench6, bench7, bench8, bench9, bench10, bench11, bench12]
N_all = [2^x for x=4:12]

using Plots, LaTeXStrings
bench_all_mean = [mean(bench_all[i].times) for i in 1:length(bench_all)]
plot(log.(2, N_all), log.(bench_all_mean), marker='*', xlab = L"n = 2^x", ylab = L"time = ln(10^9 s)", title = "log-log plaquette n x n grid", label = "naive")

fench4 = @benchmark remove_vortices_edge(findvortices_jumps(psi4), psi4)
fench5 = @benchmark remove_vortices_edge(findvortices_jumps(psi5), psi5)
fench6 = @benchmark remove_vortices_edge(findvortices_jumps(psi6), psi6)
fench7 = @benchmark remove_vortices_edge(findvortices_jumps(psi7), psi7)
fench8 = @benchmark remove_vortices_edge(findvortices_jumps(psi8), psi8)
fench9 = @benchmark remove_vortices_edge(findvortices_jumps(psi9), psi9)
fench10 = @benchmark remove_vortices_edge(findvortices_jumps(psi10), psi10)
fench11 = @benchmark remove_vortices_edge(findvortices_jumps(psi11), psi11)
fench12 = @benchmark remove_vortices_edge(findvortices_jumps(psi12), psi12)

fench_all = [fench4, fench5, fench6, fench7, fench8, fench9, fench10, fench11, fench12]
fench_all_mean = [mean(fench_all[i].times) for i in 1:length(fench_all)]

plot!(log.(2, N_all), log.(fench_all_mean), marker='*', label = "optimised")
