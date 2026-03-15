module Detection3D

abstract type AbstractDetectionBackend3D end

struct LegacyPlaneBackend3D <: AbstractDetectionBackend3D end
struct TimCopBackend3D <: AbstractDetectionBackend3D end

module TimCop

using Base.Threads
using FLoops
using Graphs
using Interpolations

using ...Analysis2D: phase_jumps

include("Detection3D/timcop.jl")

end

module Legacy

using FLoops
using Graphs
using Interpolations

using ...Analysis2D: phase_jumps

include("Detection3D/legacy.jl")

end

full_algorithm(args...; kwargs...) = TimCop.full_algorithm(args...; kwargs...)
full_algorithm(::TimCopBackend3D, args...; kwargs...) = TimCop.full_algorithm(args...; kwargs...)
full_algorithm(::LegacyPlaneBackend3D, args...; kwargs...) = Legacy.full_algorithm(args...; kwargs...)

connect_vortex_ends(args...; kwargs...) = TimCop.connect_vortex_ends(args...; kwargs...)

export AbstractDetectionBackend3D, LegacyPlaneBackend3D, TimCopBackend3D, connect_vortex_ends, full_algorithm

end
