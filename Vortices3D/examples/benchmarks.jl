using VortexDistributions, LaTeXStrings

include("../util.jl")

benchmarks = benchmark2D(200, 4, 12)

using Plots

naive = [mean(benchmarks[i][1].times) for i in 1:length(benchmarks)]
optim = [mean(benchmarks[i][2].times) for i in 1:length(benchmarks)]
optimInterp = [mean(benchmarks[i][3].times) for i in 1:length(benchmarks)]
n = [benchmarks[i][4] for i in 1:length(benchmarks)]

benchPlot = plot(log.(2, n), log.(naive), label="naive", legend=:bottomright)
plot!(log.(2, n), log.(optim), label="optimised")
plot!(log.(2, n), log.(optimInterp), label="optimisedInterpolated")
xlabel!(L"log_2(N)")
ylabel!(L"ln(10^{-9})s")
title!("Benchmark comparison of 2D Plaquette \n detection 200 vorticies")

savefig(benchPlot, "results/benchmarkPlotMultiVort200.png")

