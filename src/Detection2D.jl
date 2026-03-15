module Detection2D

using Parameters

using ..Analysis2D: phase_jumps, zoom_interp
using ..Core: Field, PointVortex, Sphere, Torus, find_where, vortex_array, Δ
using ..Creation2D: rand_vortexfield

include("Detection2D/backends.jl")
include("Detection2D/detection.jl")

export AbstractDetectionBackend2D,
    PlaquetteBackend2D,
    default_backend2d,
    findvortices,
    findvortices_grid,
    findvortices_jumps,
    found_near,
    remove_vortices_edge

end
