using VortexDistributions

## Sim params
L=(16.,16.,16.);
N=(64,64,64);
sim = Sim(L,N);
@unpack_Sim sim;