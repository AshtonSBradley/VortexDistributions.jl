##
using VortexDistributions, LaTeXStrings, Statistics, Plots
include("util.jl")
include("2Danalysis.jl")

## start at 128x128 or 256
benchmarks = benchmark2D(1, 7, 12) #128 -> 4096

##
n = [benchmarks[i][4] for i in 1:length(benchmarks)]

naive = [benchmarks[i][1] for i in 1:length(benchmarks)]
optim = [benchmarks[i][2] for i in 1:length(benchmarks)]
optimInterp = [benchmarks[i][3] for i in 1:length(benchmarks)]

mnaive = mean.(naive)
mnaive = [mnaive[i].time*10^-9 for i in 1:length(mnaive)]
σn = std.(naive)
σn = [σn[i].time*10^-9 for i in 1:length(σn)]
moptim = mean.(optim)
moptim = [moptim[i].time*10^-9 for i in 1:length(moptim)]
σo = std.(optim)
σo = [σo[i].time*10^-9 for i in 1:length(σn)]
moptimInterp = mean.(optimInterp)
moptimInterp = [moptimInterp[i].time*10^-9 for i in 1:length(moptimInterp)]
σoi = std.(optimInterp)
σoi = [σoi[i].time*10^-9 for i in 1:length(σoi)]

## 
newPlot1 = plot(log.(2, n), mnaive, ribbon=σn, yaxis=:log10)
# newPlotish = plot(log.(2, n), mnaive,  minorticks=:true, yaxis=:log10 , legend=:topleft, label="Naive", ribbon=σn)

plot!(log.(2, n), moptim, ribbon=σo, fillalpha=.4, label="Optimised")
plot!(log.(2, n), moptimInterp, ribbon=σoi, fillalpha=.4, label="Optimised w/ Interpolation")
xlabel!(L"log_2(N)")
ylabel!(L"log_{10}(t)")
title!("Benchmark comparison of 2D Plaquette \n detection 200 vorticies")
savefig(benchPlot, "results/benchmarkPlotMultiVort200.png")

## BenchmarkSuite
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 30
suite = benchmark2Dsuite(1, 7, 12)


# tune!(suite)
# BenchmarkTools.save("params.json", params(suite))
loadparams!(suite, BenchmarkTools.load("params.json")[1], :evals, :samples);
results = run(suite, verbose=true)


benchmarkDict(results)