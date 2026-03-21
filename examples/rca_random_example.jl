using Random
using VortexDistributions
using Plots

include("rca_plot_helpers.jl")

Random.seed!(31)

Lx = 2.0
Ly = 2.0
Nplus = 10
Nminus = 10

function random_vortices(n, charge; Lx = 2.0, Ly = 2.0, margin = 0.12)
    xs = rand(n) .* (Lx - 2margin) .- (Lx / 2 - margin)
    ys = rand(n) .* (Ly - 2margin) .- (Ly / 2 - margin)
    return PointVortex.(xs, ys, fill(charge, n))
end

vortices = vcat(
    random_vortices(Nplus, 1; Lx = Lx, Ly = Ly),
    random_vortices(Nminus, -1; Lx = Lx, Ly = Ly),
)

result = RCA.recursive_cluster_algorithm(vortices, Lx, Ly)
labels, colors = classify_vortices(vortices, result)

plot_path = joinpath(@__DIR__, "rca_random_example.png")
plt = plot_rca_example(
    vortices,
    result,
    labels,
    colors;
    Lx = Lx,
    Ly = Ly,
    title = "RCA Classification for Random 10+/10- Vortices with Stream Function",
)
savefig(plt, plot_path)

println("Random seed: 31")
println("Total vortices: ", length(vortices), " (", Nplus, " positive, ", Nminus, " negative)")
print_rca_summary(result, plot_path)
