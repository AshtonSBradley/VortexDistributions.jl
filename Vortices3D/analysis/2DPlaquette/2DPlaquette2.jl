## Prep
using VortexDistributions
include("Vortices3D/util.jl")
include("Vortices3D/2Danalysis.jl")

##
using Statistics
function errorSim(n, num_runs)
    Nx = n; Ny = n;
    Lx = n; Ly = Lx
    x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, Ny+1)[1:end-1];
    dx = x[2]-x[1]

    dvfn_all = []
    dvfo_all = []
    dvfoi_all = []
    for i in 1:num_runs 
        psi0 = one.(x*y') |> complex
        psi = Torus(psi0,x,y)

        rand_vort = rand_pointvortex(1, psi)
        vortex!(psi, rand_vort)
        rand_vort = rand_vort[1]

        vfnaive = naivePlaquette(psi)[1]
        vfoptim = remove_vortices_edge(findvortices_jumps(psi), psi)[1]
        vfoptimInterp = findvortices(psi)[1]

        dvfn = euclidVorts(vfnaive, rand_vort)./dx
        dvfo = euclidVorts(vfoptim, rand_vort)./dx
        dvfoi = euclidVorts(vfoptimInterp, rand_vort)./dx

        push!(dvfn_all, dvfn)
        push!(dvfo_all, dvfo)
        push!(dvfoi_all, dvfoi)
    end
    dvfn_all_mean = mean(dvfn_all)
    dvfn_all_std = std(dvfn_all)
    dvfn_all_max = findmax(dvfn_all)[1]

    dvfo_all_mean = mean(dvfo_all)
    dvfo_all_std = std(dvfo_all)
    dvfo_all_max = findmax(dvfo_all)[1]

    dvfoi_all_mean = mean(dvfoi_all)
    dvfoi_all_std = std(dvfoi_all)
    dvfoi_all_max = findmax(dvfoi_all)[1]

    return Dict([("naive", Dict([("mean", dvfn_all_mean), ("std", dvfn_all_std), ("max", dvfn_all_max)])), 
                ("vfoptim", Dict([("mean", dvfo_all_mean), ("std", dvfo_all_std), ("max", dvfo_all_max)])), 
                ("vfoptimInterp", Dict([("mean", dvfoi_all_mean), ("std", dvfoi_all_std), ("max", dvfoi_all_max)]))]);
end


result = errorSim(8, 100000) 



result["naive"]["max"]/result["vfoptimInterp"]["max"]
result["naive"]["mean"]/result["vfoptimInterp"]["mean"]

