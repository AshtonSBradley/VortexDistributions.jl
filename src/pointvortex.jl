PointVortex(v::Array{Float64,2}) = PointVortex.(v[:,1],v[:,2],v[:,3])
RawData(v::PointVortex) = [v.xv v.yv v.qv]
RawData(v::Array{PointVortex,1}) = reduce(vcat,RawData.(v))

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
function PointVortex(vort::Array{ScalarVortex{T},1}) where T
    pv = PointVortex[]
    for j in eachindex(vort)
        push!(pv,vort[j].vort)
    end
    return pv
end
