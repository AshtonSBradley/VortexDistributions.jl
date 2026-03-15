module Core

using Graphs
using Interpolations
using LinearAlgebra
using Parameters
using SimpleWeightedGraphs
using SparseArrays

include("Core/internal.jl")
include("Core/types.jl")
include("Core/pointvortex.jl")

export Ansatz,
    ChannelVortex,
    Cluster,
    CoreShape,
    Dipole,
    Exact,
    Field,
    PointVortex,
    ScalarVortex,
    Sphere,
    Torus,
    Vortex,
    VortexGroup,
    charge,
    pos,
    rand_charge,
    rand_pointvortex,
    uniform,
    vortex_array,
    xpos,
    ypos

end
