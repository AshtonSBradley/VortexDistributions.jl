"""
    PointVortex(x,y,q)

Construct a point vortex of charge `q` at `x`,`y`.

See also: [`vortex_array`](@ref), [`ScalarVortex`](@ref), [`vortex!`](@ref), [`rand_pointvortex`](@ref)
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
    varray = vortex_array(v::PointVortex)

Return `varray<:Array{Float64,2}` of vortex data.

See also: [`PointVortex`](@ref)
"""
vortex_array(v::PointVortex) = [v.xv v.yv v.qv]

function vortex_array(v::Array{PointVortex,1})
    ~isempty(v) > 0 ? (return reduce(vcat,vortex_array.(v))) : (return Array{Float64}(undef, 0, 0))
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

rand_charge() = rand([-1 1],)
rand_charge(n) = rand([-1 1],n)

"""
    rand_pointvortex()
    rand_pointvortex(n::Int)
    rand_pointvortex(ψ::Field)
    rand_pointvortex(n::Int,ψ::Field)

Sample `n::Int` random point vortices on spatial domain specified by `ψ::Field`.

See also: [`PointVortex`](@ref), [`randVortex`](@ref)
"""
rand_pointvortex() = PointVortex(uniform(),uniform(),rand_charge())
rand_pointvortex(n) = PointVortex.(uniform(n), uniform(n), rand_charge(n))

function rand_pointvortex(psi::Field)
    @unpack x,y = psi; a,b = first(x),last(x); c,d = first(y), last(y)
    return PointVortex(uniform(a,b,),uniform(c,d,),rand_charge())
end

function rand_pointvortex(n,psi::Field)
    @unpack x,y = psi;
    dx = x[2]-x[1]; dy = y[2] - y[1]
    a,b = first(x)+2dx,last(x)-2dx; c,d = first(y)+2dy, last(y)-2dy
    return PointVortex.(uniform(a,b,n),uniform(c,d,n),rand_charge(n))
end

charge(v::PointVortex) = v.qv
xpos(v::PointVortex) = v.xv
ypos(v::PointVortex) = v.yv
vpos(v::PointVortex) = [v.xv v.yv]

xpos(v::Array{PointVortex,1}) = xpos.(v)
ypos(v::Array{PointVortex,1}) = ypos.(v)
charge(v::Array{PointVortex,1}) = charge.(v)
pos(v::PointVortex)= [xpos(v) ypos(v)]
pos(v::Array{PointVortex,1}) = [xpos.(v) ypos.(v)]

pos(v::Dipole) = [xpos(v.vp) ypos(v.vp)]
