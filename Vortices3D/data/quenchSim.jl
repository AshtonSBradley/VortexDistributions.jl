using FourierGPE, VortexDistributions, GLMakie, LoopVectorization
utils = joinpath(@__DIR__, "Vortices3D/util.jl")
include(utils)
##Set simulation parameters
# L = (32., 32., 32.) #2 points per ξ
    L=(16.,16.,16.);
    N=(64,64,64);
    sim = Sim(L,N);
    @unpack_Sim sim;

## Initialse sim
    # parameters
    μ = 25.0;
    γ = 0.05;
    tf = 4/γ;
    Nt = 200;
    t = LinRange(0.,tf,Nt);

## Run sim
    x,y,z = X;
    ψi = randn(N)+im*randn(N);
    ϕi = kspace(ψi,sim);

    @pack_Sim! sim;

## Evolve in k-space
    @time sol = runsim(sim); # will take a few minutes to run.

# @save "box_vorts.jld2" psi_chaos1 psi_chaos2 psi_knots1 psi_knots2 psi_knots3 psi_knots4 psi_tubes1 psi_tubes2 psi_ringtube X

plot_iso(psi, X)



using VortexDistributions, JLD2, Interpolations;
using NearestNeighbors, LinearAlgebra, Distances;
using GLMakie, Colors
using BSON, DifferentialEquations, ScikitLearn, Static, FFTW

BSON.@save "64quench4heal.bson" sol




