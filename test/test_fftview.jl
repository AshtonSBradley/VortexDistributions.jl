#adapt FFTViews code
using Base: tail, unsafe_length, @propagate_inbounds
# using FFTW

# A custom rangetype that will be used for indices and never throws a
# boundserror because the domain is actually periodic.
using CustomUnitRanges
include(CustomUnitRanges.filename_for_urange)

@static if isdefined(Base, :IdentityUnitRange)
    const indextypes = (URange, Base.IdentityUnitRange{<:URange})
    const CircVRange{T} = Union{URange{T}, Base.Slice{URange{T}}, Base.IdentityUnitRange{URange{T}}}
    indrange(i) = Base.IdentityUnitRange(URange(first(i), last(i)))
else
    const indextypes = (URange, Base.Slice{<:URange})
    const CircVRange{T} = Union{URange{T}, Base.Slice{URange{T}}}
    indrange(i) = Base.Slice(URange(first(i), last(i)))
end

for T in indextypes
    @eval begin
        Base.checkindex(::Type{Bool}, ::$T, ::Base.Slice) = true
        Base.checkindex(::Type{Bool}, ::$T, ::Base.LogicalIndex) = true
        Base.checkindex(::Type{Bool}, ::$T, ::Real) = true
        Base.checkindex(::Type{Bool}, ::$T, ::AbstractRange) = true
        Base.checkindex(::Type{Bool}, ::$T, ::AbstractVector{Bool}) = true
        Base.checkindex(::Type{Bool}, ::$T, ::AbstractArray{Bool}) = true
        Base.checkindex(::Type{Bool}, ::$T, ::AbstractArray) = true
    end
    if isdefined(Base, :IdentityUnitRange)
        @eval Base.checkindex(::Type{Bool}, ::$T, ::Base.IdentityUnitRange) = true
    end
end

export CircView

abstract type AbstractCircView{T,N} <: AbstractArray{T,N} end

struct CircView{T,N,A<:AbstractArray} <: AbstractCircView{T,N}
    parent::A

    function CircView{T,N,A}(parent::A) where {T,N,A}
        new{T,N,A}(parent)
    end
end

CircView(parent::AbstractArray{T,N}) where {T,N} = CircView{T,N,typeof(parent)}(parent)
CircView{T,N}(dims::Dims{N}) where {T,N} = CircView(Array{T,N}(undef, dims))
CircView{T}(dims::Dims{N}) where {T,N} = CircView(Array{T,N}(undef, dims))

# Note: there are no bounds checks because it's all periodic
@inline @propagate_inbounds function Base.getindex(F::CircView{T,N}, I::Vararg{Int,N}) where {T,N}
    P = parent(F)
    @inbounds ret = P[reindex(CircView, axes(P), I)...]
    ret
end

@inline @propagate_inbounds function Base.setindex!(F::CircView{T,N}, val, I::Vararg{Int,N}) where {T,N}
    P = parent(F)
    @inbounds P[reindex(CircView, axes(P), I)...] = val
end

Base.parent(F::AbstractCircView) = F.parent
Base.axes(F::AbstractCircView) = map(indrange, axes(parent(F)))
Base.size(F::AbstractCircView) = size(parent(F))

function Base.similar(A::AbstractArray, T::Type, shape::Tuple{CircVRange,Vararg{CircVRange}})
    all(x->first(x)==1, shape) || throw(BoundsError("cannot allocate CircView with the first element of the range non-zero"))
    CircView(similar(A, T, map(length, shape)))
end

function Base.similar(f::Union{Function,Type}, shape::Tuple{CircVRange,Vararg{CircVRange}})
    all(x->first(x)==1, shape) || throw(BoundsError("cannot allocate CircView with the first element of the range non-zero"))
    CircView(similar(f, map(length, shape)))
end

Base.reshape(F::CircView{_,N}, ::Type{Val{N}}) where {_,N}   = F
Base.reshape(F::CircView{_,M}, ::Type{Val{N}}) where {_,M,N} = CircView(reshape(parent(F), Val(N)))

# FFTW.fft(F::CircView; kwargs...) = fft(parent(F); kwargs...)
# FFTW.rfft(F::CircView; kwargs...) = rfft(parent(F); kwargs...)
# FFTW.fft(F::CircView, dims; kwargs...) = fft(parent(F), dims; kwargs...)
# FFTW.rfft(F::CircView, dims; kwargs...) = rfft(parent(F), dims; kwargs...)

@inline reindex(::Type{V}, inds, I) where {V} = (_reindex(V, inds[1], I[1]), reindex(V, tail(inds), tail(I))...)
reindex(::Type{V}, ::Tuple{}, ::Tuple{}) where {V} = ()
_reindex(::Type{CircView}, ind, i) = modrange(i, ind)

modrange(i, rng::AbstractUnitRange) = mod(i-first(rng), unsafe_length(rng))+first(rng)

## test
using BenchmarkTools
x = LinRange(-5,5,1000)

@btime x[200:800]

xc = CircView(x)
@btime xc[-800:-200]


## test FFTViews separately
using FFTViews
x = LinRange(-5,5,1000)

@btime x[200:800]

xc = FFTView(x)
@btime xc[-800:-200]
