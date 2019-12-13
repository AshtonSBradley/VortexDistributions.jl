"""
    PointVortex(x,y,q)

Construct a point vortex of charge `q` at `x`,`y`.

See also: [`rawData`](@ref), [`ScalarVortex`](@ref), [`vortex!`](@ref), [`randPointVortex`](@ref)
"""
function PointVortex(v::Array{Float64,2})
    ~isempty(v) ? (return PointVortex.(v[:,1],v[:,2],v[:,3])) : (return Array{PointVortex}(undef,0))
end
function PointVortex(vort::Array{ScalarVortex{T},1}) where T
    pv = PointVortex[]
    for j in eachindex(vort)
        push!(pv,vort[j].vort)
    end
    return pv
end

"""
    varray = rawData(v::PointVortex)

Return `varray<:Array{Float64,2}` of vortex data.

See also: [`PointVortex`](@ref)
"""
rawData(v::PointVortex) = [v.xv v.yv v.qv]
function rawData(v::Array{PointVortex,1})
    ~isempty(v) > 0 ? (return reduce(vcat,rawData.(v))) : (return Array{Float64}(undef, 0, 0))
end

function uniform(a,b)
    @assert a<b
    return a + (b-a)*rand()
end
function uniform(a,b,n)
    return uniform.(a*ones(n),b)
end
uniform() = uniform(-.5,.5)
uniform(n) = uniform(-.5,.5,n)

randcharge() = rand([-1 1],)
randcharge(n) = rand([-1 1],n)

"""
    randPointVortex()
    randPointVortex(n::Int)
    randPointVortex(ψ::Field)
    randPointVortex(n::Int,ψ::Field)

Sample `n::Int` random point vortices on spatial domain specified by `ψ::Field`.

See also: [`PointVortex`](@ref), [`randVortex`](@ref)
"""
randPointVortex() = PointVortex(uniform(),uniform(),randcharge())
randPointVortex(n) = PointVortex.(uniform(n), uniform(n), randcharge(n))
function randPointVortex(psi::Field)
    @unpack x,y = psi; a,b = first(x),last(x); c,d = first(y), last(y)
    return PointVortex(uniform(a,b,),uniform(c,d,),randcharge())
end
function randPointVortex(n,psi::Field)
    @unpack x,y = psi;
    dx = x[2]-x[1]; dy = y[2] - y[1]
    a,b = first(x)+2dx,last(x)-2dx; c,d = first(y)+2dy, last(y)-2dy
    return PointVortex.(uniform(a,b,n),uniform(c,d,n),randcharge(n))
end
