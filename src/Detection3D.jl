module Detection3D

module Legacy

using FLoops
using Graphs
using Interpolations

using ...Analysis2D: phase_jumps

abstract type AbstractDetectionBackend3D end

struct LegacyPlaneBackend3D <: AbstractDetectionBackend3D end

# TODO Phase 2: replace this legacy threaded plane-linking path with the new
# VortexDetection.jl-backed 3D backend while keeping the experimental boundary.

include("Detection3D/legacy.jl")

export AbstractDetectionBackend3D, LegacyPlaneBackend3D, full_algorithm

end

end
