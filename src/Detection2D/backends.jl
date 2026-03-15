abstract type AbstractDetectionBackend2D end

struct PlaquetteBackend2D <: AbstractDetectionBackend2D end

default_backend2d(::Field) = PlaquetteBackend2D()

# TODO Phase 2: route alternative 2D detection engines through this backend boundary
# so VortexDetection.jl can plug in without reshaping the stable public API again.
