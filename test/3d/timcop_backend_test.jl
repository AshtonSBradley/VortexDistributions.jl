using Graphs
using JLD2
using Test
using VortexDistributions

@load joinpath(@__DIR__, "box_vorts.jld2") psi_tubes1 X

psi = psi_tubes1
x, y, z = X

result = VortexDistributions.detect_vortices_3d(psi, x, y, z)
g = result.graph
vort_lines = result.vortex_lines
vort_loops = result.vortex_loops
vort_rings = result.vortex_rings
vorts_coords = result.coords
edge_list_periodic = result.periodic_edges

@test Graphs.nv(g) == length(vorts_coords)
@test Graphs.ne(g) > 0
@test length(vorts_coords) > 0
@test (length(vort_lines), length(vort_loops), length(vort_rings)) == (16, 0, 1)
@test !isempty(edge_list_periodic)

connected_vorts = VortexDistributions.Detection3D.connect_vortex_ends(result, X)

@test length(connected_vorts) == 3

g2, lines2, loops2, rings2, coords2, periodic2 = VortexDistributions.full_algorithm(psi, x, y, z)
@test Graphs.nv(g2) == Graphs.nv(g)
@test Graphs.ne(g2) == Graphs.ne(g)
@test length(coords2) == length(vorts_coords)
@test length(lines2) == length(vort_lines)
@test length(loops2) == length(vort_loops)
@test length(rings2) == length(vort_rings)
@test length(periodic2) == length(edge_list_periodic)
