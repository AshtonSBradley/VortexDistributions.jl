# test typed version
# using VortexDistributions

include("typeversion.jl")

# unitary?
v0 = [.2 .4 1]
w0 = PointVortex(v0)
@test RawData(w0) == v0

v1 = [.2 .4 1;0.7 1.5 -1;-.3 1.2 1]
w1 = PointVortex(v1)
@test RawData(w1) == v1


w3 = randPointVortex(1000)
v3 = RawData(w3)
@test v3 == RawData(w3)

# fast... really, really fast!
@time v3 = RawData(w3)
@time w3 = PointVortex(v3)



function randVortexField(n)
    Nx = 400; Ny = 400
    Lx = 200; Ly = 200
    x = LinRange(-Lx / 2, Ly / 2, Nx)
    y = LinRange(-Ly / 2, Ly / 2, Ny)

    psi0 = one.(x*y') |> complex; psi = Torus(psi0,x,y)
    vort = randVortex(n,psi)
    vortex!(psi,vort)
    return psi,PointVortex(vort)
end

function PointVortex(vort::Array{ScalarVortex{T},1}) where T
    pv = PointVortex[]
    for j in eachindex(vort)
        push!(pv,vort[j].vort)
    end
    return pv
end


# Single vortex tests are reliable and don't require sorting
# (Phase cross talk will cause construction to wander)

function foundNear(n)
    near = true
    for j ∈ 1:n
        psi,vort = randVortexField(1)
        vortfound = findvortices(psi)
        vfdata = RawData(vortfound)
        vdata = RawData(vort)
        near *= isapprox(vdata,vfdata,rtol = 0.2)
    end
    return near
end

@test foundNear(1)
@test foundNear(10)
@test foundNear(30)

psi,vort = randVortexField(1);heatmap(psi.x,psi.y,abs2.(psi.ψ))

# test real data with masking examples

@load "./test/one_frame.jld2"
x0,y0 = 70,200

# circular
msk(x,y) = sqrt(x^2 + y^2) > 100
@. ψ1[msk(x-x0,y'-y0) ] = false
heatmap(x,y,abs2.(ψ1),transpose=true)

# square
Lx,Ly = 100,100

xm = x[@. abs(x-x0)<Lx/2]
ym = y[@. abs(y-y0)<Ly/2]
psim = @. ψ1[abs(x-x0)<Lx/2,abs(y-y0)<Ly/2]

heatmap(xm,ym,angle.(psim),transpose=true)
psi2 = Torus(psim,xm,ym)
vort = findvortices(psi2)

scatter!([vort[1].xv],[vort[1].yv],alpha=0.4)
scatter!([vort[2].xv],[vort[2].yv],alpha=0.4)
