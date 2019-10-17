# make detection more robust
# 2. sub-grid detection based on fitting, rather than accurate zero finding
# rationale: core location isn't a good proxy for vortex location.
# vortex is a bulk flow excitation, whereas any small noise will shift zero - without shifting bulk

#---
using Plots, Revise, VortexDistributions

# unitary?
v0 = [.2 .4 1]
w0 = PointVortex(v0)
@test rawData(w0) == v0

# @test foundNear(1)

# no vortices
Nx = 400; Ny = 400
Lx = 200;
Ly = 200;
x = LinRange(-Lx / 2, Ly / 2, Nx)
y = LinRange(-Ly / 2, Ly / 2, Ny)

psi0 = one.(x * y') |> complex;
psi = Torus(psi0, x, y);
vortfound = findvortices(psi)
@test vortfound == nothing


#--- debugging functions
import VortexDistributions: findvortices, findvortices_grid, findvortices_interp, findvortices_jumps

function findvortices(psi::Field,interp=true)
    if interp
        return findvortices_interp(psi)
    else
        return findvortices_grid(psi)
    end
end

function findvortices_grid(psi::Torus;shift=true)
    vort = findvortices_jumps(psi,shift=shift)
    rawvort = rawData(vort)
    vort = sortslices(rawvort,dims=1)
    return PointVortex(vort)
end

function findvortices_jumps(psi::Field;shift=true)
    @unpack x,y,ψ = psi
    phase = angle.(ψ)

    # Compute all nearest neighbour jumps
    diffx = phasejumps(phase,1); diffy = phasejumps(phase,2)

    # Compute all plaquette loop integrals with in-place memory recycling
    circshift!(phase,diffx,(0,1))
    diffx .-= phase; diffx .-= diffy
    circshift!(phase,diffy,(1,0))
    diffx .+= phase

    # Find all windings
    ixp,iyp,vp = findwhere(diffx .> 0.0)
    xp = x[ixp]; yp = y[iyp]; np = length(vp)
    ixn,iyn,vn = findwhere(diffx .< 0.0)
    xn = x[ixn]; yn = y[iyn]; nn = length(vn)

    if shift
        dx = x[2]-x[1]; dy = y[2] - y[1]
        xp .-= dx/2; yp .-= dy/2; xn .-= dx/2; yn .-= dy/2
    end

    vortices = [xn yn -vn; xp yp vp]

    return PointVortex(vortices)
end

function findvortices_interp(psi::Field)
    vort = findvortices_grid(psi,shift=true)

        return vort
    end
end

# function findvortices_interp(psi::Field)
#     vort = findvortices_grid(psi,shift=true)
#     vort = removeedgevortices(vort,psi)
#     #TODO: allow for interp with periodic data (here edges are stripped)
#
#     for (j,vortex) in enumerate(vort)
#         try
#         vortz,psiz = corezoom(vortex,psi)
#         vortz,psiz = corezoom(vortz,psiz)
#         vortz,psiz = corezoom(vortz,psiz)
#         vort[j] = vortz
#         catch nothing
#         end
#         return vort
#     end
# end

#--- one vortex
Nx = 400; Ny = 400
Lx = 200;
Ly = 200;
x = LinRange(-Lx / 2, Ly / 2, Nx)
y = LinRange(-Ly / 2, Ly / 2, Ny)

psi0 = one.(x * y') |> complex;
psi = Torus(psi0, x, y);
vort = randVortex(1, psi)
vortex!(psi, vort)

vort = PointVortex(vort)


vortfound = findvortices(psi)

vfdata = rawData(vortfound)

vdata = rawData(vort)

heatmap(x,y,abs2.(psi.ψ),aspect_ratio=1)
heatmap(x,y,angle.(psi.ψ),aspect_ratio=1)


# Let's try fitting a Gaussian to the density and compare with root finding
# Do this with and without local grid noise

#--- local field (within 2 healing lengths)
x1 = x[@. abs(x -vort[1].xv)<2]
y1 = y[@. abs(x -vort[1].yv)<2]

psi1 = psi.ψ[(@. abs(x -vort[1].xv)<2) , @. abs(x -vort[1].yv)<2]

heatmap(x1,y1,abs2.(psi1),aspect_ratio=1)
heatmap(x1,y1,angle.(psi1),aspect_ratio=1)

m = abs2.(psi1)
m = maximum(m) .- m
heatmap(x1,y1,m,aspect_ratio=1)

using DataFitting

gaussian(x, y, p1, p2, p3,p4) = @.  exp(-((x-p1)^2  +  (y-p2)^2)/p3^2)*p4

dom = CartesianDomain(x1,y1)

data = Measures(m) #Measures(d + randn(rng, size(d)), 1.)

model1 = Model(:comp1 => FuncWrap(gaussian, vort[1].xv, vort[1].yv, 1.,1.))
prepare!(model1, dom, :comp1)
@time result = fit!(model1, data)


#--- add local noise
phi = Torus(copy(psi.ψ),x,y)
phi.ψ .+= (randn()+im*randn())*0.2

vort = findvortices_interp(phi)
x1 = x[@. abs(x -vort[1].xv)<2]
y1 = y[@. abs(x -vort[1].yv)<2]

phi1 = phi.ψ[(@. abs(x -vort[1].xv)<2) , @. abs(x -vort[1].yv)<2]

heatmap(x1,y1,abs2.(phi1),aspect_ratio=1)
heatmap(x1,y1,angle.(phi1),aspect_ratio=1)
heatmap(x,y,angle.(phi.ψ),aspect_ratio=1)

m = abs2.(phi1)
m = maximum(m) .- m
heatmap(x1,y1,m,aspect_ratio=1)

using DataFitting

gaussian(x, y, p1, p2, p3,p4) = @.  exp(-((x-p1)^2  +  (y-p2)^2)/p3^2)*p4

dom = CartesianDomain(x1,y1)

data = Measures(m) #Measures(d + randn(rng, size(d)), 1.)

model1 = Model(:comp1 => FuncWrap(gaussian, vort[1].xv, vort[1].yv, 1.,1.))
prepare!(model1, dom, :comp1)
@time result = fit!(model1, data)

vfoundnoise = findvortices(phi) |> rawData
