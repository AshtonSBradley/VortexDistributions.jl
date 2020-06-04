abstract type Field end
abstract type Vortex end
abstract type CoreShape end
abstract type VortexGroup end

"""
    Torus <: Field

Specify toroidal boundary conditions: doubly periodic in 2D.

# Arguments
- `ψ::Array{Complex{Float64},2}`: wavefunction.
- `x::Vector{Float64}`: `x` coordinates.
- `y::Vector{Float64}`: `y` coordinates.

See also: [`Sphere`](@def), [`Field`](@ref)
"""
mutable struct Torus <: Field
    ψ::Array{Complex{Float64},2}
    x::Vector{Float64}
    y::Vector{Float64}
end

"""
    Sphere <: Field

Specify spherical boundary conditions in 2D.

# Arguments
- `ψ::Array{Complex{Float64},2}`: wavefunction.
- `x::Vector{Float64}`: `x` coordinates.
- `y::Vector{Float64}`: `y` coordinates.

See also: [`Torus`](@def), [`Field`](@ref)
"""
mutable struct Sphere <: Field
    ψ::Array{Complex{Float64},2}
    x::Vector{Float64}
    y::Vector{Float64}
end

mutable struct PointVortex <: Vortex
    x::Float64
    y::Float64
    q::Int64
end

struct Ansatz <: CoreShape
    f::Interpolations.GriddedInterpolation{Float64,1,Float64,Gridded{Linear},Tuple{Array{Float64,1}}}
    ξ::Float64
    Λ::Float64
end

struct Exact <: CoreShape
    f::Interpolations.GriddedInterpolation{Float64,1,Float64,Gridded{Linear},Tuple{Array{Float64,1}}}
    ξ::Float64
end

mutable struct ScalarVortex{T<:CoreShape} <: Vortex
    core::T
    vort::PointVortex
end

mutable struct Cluster <: VortexGroup
    v::Array{PointVortex,1}
end

mutable struct Dipole <: VortexGroup
    vp::PointVortex
    vn::PointVortex
    Dipole(v1,v2) = (charge(v1) > 0 && charge(v2) < 0) ? new(v1,v2) : new(v2,v1)
end
