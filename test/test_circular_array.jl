struct CircularArray{T,N} <: AbstractArray{T,N} #inherits from AbstractArray
    x::AbstractArray{T,N}
    function CircularArray(x::AbstractArray{T,N}) where {T,N}
#creates the type only with a vector x
        return new{T,N}(x)
    end
end

Base.size(A::CircularArray) = size(A.x) #important: use of Base.function

Base.length(A::CircularArray)=length(A.x)

function Base.getindex(A::CircularArray, I::Vararg{Int, N}) where N # implements A[I]
    I2 = size(A)
    return Base.getindex(A.x,(mod.(I .- 1,I2) .+ 1)...) #this is the magic operation
end

Base.getindex(A::CircularArray, I) = (A[i] for i in I) #A[1:5], for example

function Base.setindex!(A,value,I::Vararg{Int, N}) where N # A[I] = value
    return Base.setindex!(A.x,value,(mod.(I .- 1,I2) .+ 1)...)
end

# a fix for display issue?
function Base.setindex!(::Dict{Int64,V},::Any,::Int64) where V

Base.IndexStyle(::Type{CircularArray}) = IndexCartesian()

Base.checkbounds(A::CircularArray, I...) = nothing

## test
using BenchmarkTools
x = LinRange(-5,5,1000)

@btime x[200:800]

xc = CircularArray(x)
xnew = xc[-3:-1]
