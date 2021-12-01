using VortexDistributions 

Nx=1000; Ny=1000;
θx = LinRange(0,4*pi,Nx); θy = LinRange(0,4*pi,Ny);
psix = exp.(im*θx) + 0.1*(randn(Nx)+im*randn(Nx)); psiy = exp.(im*θy) + 0.1*(randn(Ny)+im*randn(Ny));
ϕx = angle.(psix); ϕy = angle.(psiy);
ϕ = [ϕx ϕy]
ϕu = unwrap(ϕ)
ϕui = similar(ϕu)
unwrap!(ϕui,ϕ)
@test ϕui == ϕu

ϕ = zeros(2, Nx)
ϕ[1, :] .= ϕx
ϕ[2, :] .=ϕy
ϕu = unwrap(ϕ,2)
ϕui = similar(ϕu)
unwrap!(ϕui,ϕ,2)
@test ϕui == ϕu
