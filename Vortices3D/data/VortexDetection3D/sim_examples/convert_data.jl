using FourierGPE
using VortexDistributions
using JLD2

@load "src/VortexDetection3D/sim_examples/3dquenchslab_data.jld2" psi1 psi2 psi3

## parameters
L = (30.,30.,8.)
N = (128,128,32)
sim = Sim(L,N)
@unpack_Sim sim;


## Initialize simulation
γ = 0.05
μ = 25.0
tf = 2.5/γ
Nt = 100
t = LinRange(0.,tf,Nt)


## random initial state
x,y,z = X
ψi = randn(N)+im*randn(N);
ϕi = kspace(ψi,sim);

@pack_Sim! sim;

psi_tubes1 = xspace(psi1, sim)
psi_tubes2 = xspace(psi2, sim)
psi_tubes3 = xspace(psi3, sim)

@save "src/VortexDetection3D/sim_examples/slab_tubes.jld2" psi_tubes1 psi_tubes2 psi_tubes3 X