using FourierGPE
using VortexDistributions
using GraphMakie
using GLMakie


## Load sim
@load "sim64.jld2" sim
@load "sol64.jld2" sol

## params
t_idx = 27
psi = sol[t_idx]
X = sim.X;
x = X[1]; y = X[2]; z = X[3];
x = round.(x, digits=3); y = round.(y, digits=3); z = round.(z, digits=3);
dx = x[2] - x[1]; dy = y[2] - y[1]; dz = z[2] - z[1];

## Isosurface
density = abs2.(psi)
pmax = maximum(density)
density = density/pmax
xlims = (X[1][1], X[1][end])
ylims = (X[1][1], X[1][end])
zlims = (X[1][1], X[1][end])
volume(xlims, ylims, zlims, density, algorithm = :iso, transparency=true, alpha=0.2)

## Detection


@time g, vort_lines, vort_loops, vort_rings, vorts_coords = full_algorithm(psi, x, y, z, n_itp = 4);
graphplot!(g, layout = (adj) -> vorts_coords, markersize=0.02, edge_width=1, node_size=5)

## Smoothing

vort_lines_coords = [reduce(vcat, transpose.(vorts_coords[v_line])) for v_line in vort_lines];
vort_loops_coords = [reduce(vcat, transpose.(vorts_coords[v_loop])) for v_loop in vort_loops];
vort_rings_coords = [reduce(vcat, transpose.(vorts_coords[v_ring])) for v_ring in vort_rings];

@time vort_lines_ccma_j, vort_loops_ccma_j, vort_rings_ccma_j = vortex_ccma_j(vort_lines_coords, vort_loops_coords, vort_rings_coords, w_ma = 6, w_cc = 3, distrib = "hanning");

# graphplot(g, layout = (adj) -> vorts_coords, markersize=0.02, edge_width=0, node_size=5)
plot_unconnected_vorts_mat(vort_lines_ccma_j, vort_loops_ccma_j, vort_rings_ccma_j, linewidth=5)