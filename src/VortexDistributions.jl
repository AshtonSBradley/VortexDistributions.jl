module VortexDistributions

using Test
using JLD2
using Parameters
using SpecialFunctions
using Interpolations
using LinearAlgebra
using ToeplitzMatrices
using SparseArrays
using FFTW
using FileIO
using ProgressMeter
using LightGraphs
using SimpleWeightedGraphs

const Λ = 0.8249

# topology
export Field, Torus, Sphere,

# cores
Vortex, CoreShape, Ansatz, Exact, ScalarVortex, PointVortex,

# vortex groups
Dipole, Cluster, VortexGroup, 

# detection 
findvortices, found_near, phase_jumps, phase_jumps!, unwrap, unwrap!, Δ,
find_where, findvortices_jumps, findvortices_grid,
remove_vortices_edge, zoom_interp, zoom_grid,

# construction
scalar_ansatz, vortex_array, uniform,
vortex!, dipole_phase, periodic_dipole!,
rand_charge, rand_pointvortex, rand_scalarvortex, rand_vortex, rand_vortexfield,

# convenient access
charge, xpos, ypos, pos

include("RecursiveClusterAlgothm/RecursiveClusterAlgorithm.jl")

include("types.jl")
include("pointvortex.jl")
include("detection.jl")
include("creation.jl")
include("utils.jl")

@load joinpath(@__DIR__,"exactcore.jld2") ψi
@load joinpath(@__DIR__,"ansatzcore.jld2") ψa

end
