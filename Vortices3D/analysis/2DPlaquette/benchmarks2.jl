##
using VortexDistributions, Statistics
include("Vortices3D/util.jl")
include("Vortices3D/2Danalysis.jl")

## Params
Nstart = 7; Nend = 12 # 2^7 = 128 -> 2^12 = 4096
num_vorts = 60; # Psi with one vortex for each trial
maxSec = 10 # Default is 5sec which is too small for larger grids which take ~2sec

## Run benchmark 
suite = benchmark2Dsuite(maxSec, num_vorts, Nstart, Nend) # Creates the test suite 
result60 = run(suite, verbose=true) # Runs the suite
results60 = benchmarkDict(result)  # Takes the output dict and formats it for analysis
results1 = benchmarkDict(result)
results10 = benchmarkDict(result10)
results30 = benchmarkDict(result30)
results60 = benchmarkDict(result60)

using JLD2
@save "benchSuite10vort360sec.jld2" result

test = copy(result)
@load "benchSuite10vort360sec.jld2" result

## Plotting
using Plots, LaTeXStrings
plots = []
algo_keys = ["naive", "optim", "optimInterp"] 

x = results1["naive"]["min"][:, 1]; 
y1 = results1["naive"]["min"][:, 2];
y2 = results1["optim"]["min"][:, 2];
y3 = results1["optimInterp"]["min"][:, 2];
p = plot(log.(2, x), y1*10^-9, label="naive", yscale=:log10, legend=:topleft, minorticks=:true)
plot!(log.(2, x), y2*10^-9, label="optim")
plot!(log.(2, x), y3*10^-9, label="optimInterp")
xlabel!(L"log_2(N)")
ylabel!(L"log_{10}(t) \enspace [s]")
push!(plots, p)

y1 = results10["naive"]["min"][:, 2];
y2 = results10["optim"]["min"][:, 2];
y3 = results10["optimInterp"]["min"][:, 2];
p = plot(log.(2, x), y1*10^-9, label="naive", yscale=:log10, legend=:topleft, minorticks=:true)
plot!(log.(2, x), y2*10^-9, label="optim")
plot!(log.(2, x), y3*10^-9, label="optimInterp")
xlabel!(L"log_2(N)")
ylabel!(L"log_{10}(t) \enspace [s]")
push!(plots, p)

y1 = results30["naive"]["min"][:, 2];
y2 = results30["optim"]["min"][:, 2];
y3 = results30["optimInterp"]["min"][:, 2];
p = plot(log.(2, x), y1*10^-9, label="naive", yscale=:log10, legend=:topleft, minorticks=:true)
plot!(log.(2, x), y2*10^-9, label="optim")
plot!(log.(2, x), y3*10^-9, label="optimInterp")
xlabel!(L"log_2(N)")
ylabel!(L"log_{10}(t) \enspace [s]")
push!(plots, p)

y1 = results60["naive"]["min"][:, 2];
y2 = results60["optim"]["min"][:, 2];
y3 = results60["optimInterp"]["min"][:, 2];
p = plot(log.(2, x), y1*10^-9, label="naive", yscale=:log10, legend=:topleft, minorticks=:true)
plot!(log.(2, x), y2*10^-9, label="optim")
plot!(log.(2, x), y3*10^-9, label="optimInterp")
xlabel!(L"log_2(N)")
ylabel!(L"log_{10}(t) \enspace [s]")
push!(plots, p)

combined = plot(plots...)
savefig(combined, "Vortices3D/results/benchAllVorts.png")

## Ratio plots

# x = results1["naive"]["min"][:, 1]; 
# y1 = results1["naive"]["min"][:, 2];
# y2 = results1["optim"]["min"][:, 2];
# y3 = results1["optimInterp"]["min"][:, 2];

# plot(log.(2, x), ones(length(x)), label="naive", legend=:topleft, yscale=:log, minorticks=:true)
# plot!(log.(2, x), y2./y1, label="optim")
# plot!(log.(2, x), y3./y1, label="optimInterp")



