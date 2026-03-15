function PointVortex(v::Array{Float64, 2})
    return isempty(v) ? PointVortex[] : PointVortex.(v[:, 1], v[:, 2], v[:, 3])
end

function PointVortex(vort::Array{ScalarVortex{T}, 1}) where {T}
    pv = PointVortex[]
    for j in eachindex(vort)
        push!(pv, vort[j].vort)
    end
    return pv
end

vortex_array(v::PointVortex) = [v.xv v.yv v.qv]

function vortex_array(v::Array{PointVortex, 1})
    return isempty(v) ? Array{Float64}(undef, 0, 0) : reduce(vcat, vortex_array.(v))
end

function uniform(a, b)
    @assert a < b
    return a + (b - a) * rand()
end

uniform(a, b, n) = uniform.(a * ones(n), b)
uniform() = uniform(-0.5, 0.5)
uniform(n) = uniform(-0.5, 0.5, n)

rand_charge() = rand([-1 1],)
rand_charge(n) = rand([-1 1], n)

rand_pointvortex() = PointVortex(uniform(), uniform(), rand_charge())
rand_pointvortex(n) = PointVortex.(uniform(n), uniform(n), rand_charge(n))

function rand_pointvortex(psi::Field)
    @unpack x, y = psi
    a, b = first(x), last(x)
    c, d = first(y), last(y)
    return PointVortex(uniform(a, b), uniform(c, d), rand_charge())
end

function rand_pointvortex(n, psi::Field)
    @unpack x, y = psi
    dx = x[2] - x[1]
    dy = y[2] - y[1]
    a, b = first(x) + 2dx, last(x) - 2dx
    c, d = first(y) + 2dy, last(y) - 2dy
    return PointVortex.(uniform(a, b, n), uniform(c, d, n), rand_charge(n))
end

charge(v::PointVortex) = v.qv
xpos(v::PointVortex) = v.xv
ypos(v::PointVortex) = v.yv
vpos(v::PointVortex) = [v.xv v.yv]

xpos(v::Array{PointVortex, 1}) = xpos.(v)
ypos(v::Array{PointVortex, 1}) = ypos.(v)
charge(v::Array{PointVortex, 1}) = charge.(v)
pos(v::PointVortex) = [xpos(v) ypos(v)]
pos(v::Array{PointVortex, 1}) = [xpos.(v) ypos.(v)]

pos(v::Dipole) = [xpos(v.vp) ypos(v.vp)]
