##
using Plots, VortexDistributions

##
N=50
θ = LinRange(0,2*pi,N)
psi = exp.(im*θ) 
ϕ = angle.(psi)
@time ϕu = unwrap(ϕ)
ϕui = similar(ϕu)

@time unwrap!(ϕui,ϕ)

## 
plot(ϕ,label="phi",c=:red,alpha=0.4,lw=6)
plot!(ϕu,label="phiu",lw=2,c=:black)
plot!(ϕui,label="phiu",lw=2,c=:green)