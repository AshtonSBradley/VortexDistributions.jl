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

function rawdata(v::Array{PointVortex,1})
    return [getx.(v) gety.(v) getq.(v)]
end

@time varray = randvortex(10000)
@time vd = rawdata(varray)
