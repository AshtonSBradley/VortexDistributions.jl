using VortexDistributions, JLD2, Interpolations;
using NearestNeighbors, LinearAlgebra, Distances;
using GLMakie, Colors

include("../util.jl");
# @load "Vortices3D/data/box_vorts.jld2" psi_chaos1 psi_chaos2 psi_knots1 psi_knots2 psi_knots3 psi_knots4 psi_tubes1 psi_tubes2 psi_ringtube X
@load "data/box_rings.jld2" psi_ring1 psi_ring2 psi_ring3 X # make sure julia repl is in base directory 


psi = psi_ring1;
plot_iso(psi,X)

## Find vortex points
N = 100
@time vorts_3d = findvortices3D_itp(psi, X, N); # Find vortex points with interpolation depth N
v_matrix = vcat(vorts_3d'...)[:,1:3]' # Convert to matrix for kdtree 

scatterVortsOnIso(vorts_3d)
## Find filaments
x = X[1]; y = X[2]; z = X[3];
Δx = x[2]-x[1]; Δy = y[2]-y[1]; Δz = z[2]-z[1];

α = .2
ϵ = (1+α)*sqrt(Δx^2+Δy^2+Δz^2)/N
@time sets3 = setMethod3(v_matrix, ϵ)
@time setsPeriodic = setMethodPeriodic(v_matrix, X, ϵ, true)

# @time nearest = uniqueVortices(ϵ, psi, X, N)
# @time sets1 = setMethod(psi, X, N, ϵ)
# @time sets2 = setMethod2(v_matrix, ϵ);


plot_iso(psi, X)
using Colors
colors = distinguishable_colors(length(sets3),[RGB(0,0.35,0.25)],dropseed=true)
for i in 1:length(setsPeriodic)
    vi = v_matrix[:, collect(setsPeriodic[i])]
    scatter!(vi[1,:],vi[2,:],vi[3,:],markersize=200,color=colors[i])
end

sortslices(v_matrix, dims=2, by=x->x[1])

A = [v_matrix[:, i] for i in 1:length(v_matrix[1,:]) if v_matrix[1, i] == -7.75 || v_matrix[1,i] == 8.0]

for i in 1:length(A)
    println(A[i])
end

## x[end] points
scatter!([8.0], [-0.875], [1.2947], markersize=400, color="red") #A[end-9]

scatter!([8.0], [-3.34913], [0.741825], markersize=400, color="red") #A[end-11]

## x[1] points 

scatter!([-7.75], [-1.03938], [0.9677467], markersize=400, color="blue") #A[9]

scatter!([-7.75], [-3.707936], [0.317925], markersize=400, color="blue") #A[6]

euclidean(A[end-9][2:3], A[9][2:3])
euclidean(A[end-11][2:3], A[6][2:3])

sqrt(3)Δx



@time psi_itp1 = interpolate(psi, BSpline(Cubic(Periodic(OnCell()))));
@time psi_itp2 = interpolate(psi, BSpline(Quadratic(Periodic(OnCell()))));


psi_exp = extrapolate(psi_itp, Periodic(OnCell()))

psi_exp[0.1,0:65,0:65]

psi_itp[1.1,1.1,1.1]
psi_exp[65.1,65.1,65.1]

x_itp = interpolate(X[1], BSpline(Linear()))
x_etp = extrapolate(x_itp, Line())
x_etp[64]

x_etp[1:65]

r1 = LinRange(1,10,10)
r2 = LinRange(10,11,10)
vcat(r1, r2)

X[1]

x = collect(LinRange(1,2,64))
X[1]