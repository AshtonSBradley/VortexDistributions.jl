using VortexDistributions, Plots

N=1000
θ = LinRange(0,4*pi,N)
psi = exp.(im*θ) + 0.1*(randn(N)+im*randn(N))
ϕ = angle.(psi)
ϕu = unwrap(ϕ)
ϕui = similar(ϕu)
unwrap!(ϕui,ϕ)

@test ϕu == ϕui


# plot(ϕ,label="phi",c=:red,alpha=0.4,lw=6)
# plot!(ϕu,label="phiu",lw=2,c=:black)
# plot!(ϕui,label="phiu",lw=2,c=:green)