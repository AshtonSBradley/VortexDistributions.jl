abstract type Field end
abstract type Vortex end
abstract type CoreShape end

@with_kw mutable struct Torus <: Field
    ψ::Array{Complex{Float64},2}
    x::Vector{Float64}
    y::Vector{Float64}
end

@with_kw mutable struct Sphere <: Field
    ψ::Array{Complex{Float64},2}
    x::Vector{Float64}
    y::Vector{Float64}
end

@with_kw mutable struct PointVortex <: Vortex
    xv::Float64
    yv::Float64
    qv::Int64
end

@with_kw mutable struct Ansatz <: CoreShape
    f::Interpolations.GriddedInterpolation{Float64,1,Float64,Gridded{Linear},Tuple{Array{Float64,1}}}
    ξ::Float64
    Λ::Float64
end
@with_kw mutable struct Exact <: CoreShape
    f::Interpolations.GriddedInterpolation{Float64,1,Float64,Gridded{Linear},Tuple{Array{Float64,1}}}
    ξ::Float64
end
@with_kw mutable struct ScalarVortex{T<:CoreShape} <: Vortex
    core::T
    vort::PointVortex
end
