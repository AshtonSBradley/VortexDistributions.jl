module Creation2D

using FFTW
using Interpolations
using JLD2
using LinearAlgebra
using Parameters
using SpecialFunctions

using ..Core: CoreShape, Field, PointVortex, Torus, Vortex, rand_pointvortex, vortex_array
import ..Core: Ansatz, ChannelVortex, Exact, ScalarVortex

const Λ = 0.8249

JLD2.@load joinpath(@__DIR__, "cores.jld2") ψi ψa

include("Creation2D/creation.jl")

export dipole_phase,
    periodic_dipole!,
    rand_charge,
    rand_scalarvortex,
    rand_vortex,
    rand_vortexfield,
    scalar_ansatz,
    thetad,
    vortex!

using ..Core: rand_charge

end
