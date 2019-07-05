# make and save charge 1 core for later interp
using Plots, LaTeXStrings, VortexDistributions

# n=1
# y,ψ,res = gpecore(n)
# @save "./src/vortexcore.jld2" y ψ
#
#
# loadpath = dirname(pathof(VortexDistributions))
# @load loadpath*"/vortexcore.jld2" y ψ
#     ψi = interpolate(tuple(y[1:end-1]), ψ[1:end-1], Gridded(Linear()))

n(x) = x^2/(1+x^2)
v(x) = 1/x
j(x) = n(x)*v(x)
nv(x) = vortexcore(x,1,false)^2
jv(x) = nv(x)*v(x)

nx = 1000
xmax = 6
x = LinRange(0,xmax,nx)
plot(x,nv.(x),label=L"n_v(r)",fill=(0,0.3))
plot!(x,n.(x),label=L"n_a(r)")
xlabel!(L"r/\xi")
ylabel!(L"n/n_0")
