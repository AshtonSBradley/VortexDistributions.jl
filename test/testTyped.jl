# test typed version
include("typeversion.jl")
v1 = PointVortex(.1,.2,1)
v2 = PointVortex(randn(),randn(),rand([1,-1]))

function randvortex(n)
    return PointVortex.(randn(n),randn(n),rand([1,-1],n))
end

getx(v::PointVortex) = (p->p.xv)(v)
gety(v::PointVortex) = (p->p.yv)(v)
getq(v::PointVortex) = (p->p.qv)(v)

function rawvortices(v::Array{PointVortex,1})
    return [getx.(v) gety.(v) getq.(v)]
end

@time varray = randvortex(10000)
@time vd = rawvortices(varray)

function PointVortex(v::Array{Float64,2})
    @assert size(v) == (1,3)
    return PointVortex(v...)
end

function pointvortices(v::Array{Float64,2})
    @assert size(v)[2] == 3
    vort = Array{PointVortex,1}(undef,0)
    for j in 1:size(v)[1]
        push!(vort,PointVortex(v[j,:]))
    end
    return vort
end

v0 = [.2 .4 1]

w0 = PointVortex(v0)

v1 = [.2 .4 1;0.7 1.5 -1;-.3 1.2 1]

w1 = pointvortices(v1)

v2 = rawvortices(w1)

v2 == v1
