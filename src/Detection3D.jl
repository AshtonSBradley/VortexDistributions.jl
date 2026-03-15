module Detection3D

abstract type AbstractDetectionBackend3D end

struct LegacyPlaneBackend3D <: AbstractDetectionBackend3D end
struct TimCopBackend3D <: AbstractDetectionBackend3D end

struct Detection3DResult{G, C, E}
    graph::G
    vortex_lines::Vector
    vortex_loops::Vector
    vortex_rings::Vector
    coords::C
    periodic_edges::E
end

function Detection3DResult(graph, vortex_lines, vortex_loops, vortex_rings, coords;
                           periodic_edges = Tuple{Int, Int}[])
    return Detection3DResult(graph, vortex_lines, vortex_loops, vortex_rings, coords, periodic_edges)
end

module TimCop

using Base.Threads
using FLoops
using Graphs
using Interpolations

using ...Analysis2D: phase_jumps
using ..Detection3D: Detection3DResult

include("Detection3D/timcop.jl")

end

module Legacy

using FLoops
using Graphs
using Interpolations

using ...Analysis2D: phase_jumps

include("Detection3D/legacy.jl")

end

detect_vortices_3d(args...; backend::AbstractDetectionBackend3D = TimCopBackend3D(), kwargs...) =
    detect_vortices_3d(backend, args...; kwargs...)
detect_vortices_3d(::TimCopBackend3D, args...; kwargs...) = TimCop.detect_vortices_3d(args...; kwargs...)
detect_vortices_3d(::LegacyPlaneBackend3D, args...; kwargs...) = Detection3DResult(Legacy.full_algorithm(args...; kwargs...)...)

function full_algorithm(args...; backend::AbstractDetectionBackend3D = TimCopBackend3D(), kwargs...)
    result = detect_vortices_3d(args...; backend = backend, kwargs...)
    return result.graph,
        result.vortex_lines,
        result.vortex_loops,
        result.vortex_rings,
        result.coords,
        result.periodic_edges
end

full_algorithm(::TimCopBackend3D, args...; kwargs...) = full_algorithm(args...; backend = TimCopBackend3D(), kwargs...)
full_algorithm(::LegacyPlaneBackend3D, args...; kwargs...) = full_algorithm(args...; backend = LegacyPlaneBackend3D(), kwargs...)

connect_vortex_ends(args...; kwargs...) = TimCop.connect_vortex_ends(args...; kwargs...)
connect_vortex_ends(result::Detection3DResult, X) =
    connect_vortex_ends(result.coords, result.vortex_lines, result.vortex_loops, result.vortex_rings, X)

export AbstractDetectionBackend3D,
    Detection3DResult,
    LegacyPlaneBackend3D,
    TimCopBackend3D,
    connect_vortex_ends,
    detect_vortices_3d,
    full_algorithm

end
