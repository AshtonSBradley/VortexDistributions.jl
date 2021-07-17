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
    xv::Float64
    yv::Float64
    qv::Int64
end

struct Ansatz <: CoreShape
    f::Interpolations.GriddedInterpolation{Float64, 1, Float64, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float64}}}
    ξ::Float64
    Λ::Float64
end

struct Exact <: CoreShape
    f::Interpolations.GriddedInterpolation{Float64, 1, Float64, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float64}}}
    ξ::Float64
end

mutable struct ScalarVortex{T<:CoreShape} <: Vortex
    core::T
    vort::PointVortex
end

mutable struct ChannelVortex{T<:CoreShape} <: Vortex
    core::T
    vort::PointVortex
end


mutable struct Dipole <: VortexGroup
    vp::PointVortex
    vn::PointVortex
    Dipole(v1,v2) = (charge(v1) > 0 && charge(v2) < 0) ? new(v1,v2) : new(v2,v1)
end

mutable struct Cluster <: VortexGroup
    vortices::Array{PointVortex,1}
    tree::Array{LightGraphs.SimpleGraphs.SimpleEdge{Int64},1}
end
Cluster(vort::Array{PointVortex,1}) = Cluster(vort,spanning_tree(xpos(vort),ypos(vort)))
