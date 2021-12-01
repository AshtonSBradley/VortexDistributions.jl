using VortexDistributions 

N=1000
θ = LinRange(0,4*pi,N)
psi = exp.(im*θ) + 0.1*(randn(N)+im*randn(N))
ϕ = angle.(psi)
ϕu = unwrap(ϕ)
ϕui = similar(ϕu)
unwrap!(ϕui,ϕ)

@test ϕu == ϕui
