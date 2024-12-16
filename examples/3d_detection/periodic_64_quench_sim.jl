using FourierGPE
using VortexDistributions

## Creating a 3D simulation

L=(16.,16.,16.);
N=(64,64,64);
sim = Sim(L,N);
@unpack_Sim sim;

μ = 25.0;
γ = 0.05;
tf = 4/γ;
Nt = 50;
t = LinRange(0.,tf,Nt);

## Run sim
x,y,z = X;
ψi = randn(N)+im*randn(N);
ϕi = kspace(ψi,sim);

@pack_Sim! sim;

## Evolve in k-space
# import FourierGPE
@time sol = runsim(sim); # will take a few minutes to run.

# Mutate the solution to xspace before saving to jld2
for i in eachindex(sol)
    sol[i] = xspace(sol[i], sim)
end

@save "sol64.jld2" sol
@save "sim64.jld2" sim