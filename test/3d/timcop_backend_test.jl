using Graphs
using JLD2
using Test
using VortexDistributions

@load joinpath(@__DIR__, "box_vorts.jld2") psi_tubes1 X

psi = psi_tubes1
x, y, z = X

g, vort_lines, vort_loops, vort_rings, vorts_coords, edge_list_periodic =
    VortexDistributions.Detection3D.TimCop.full_algorithm(psi, x, y, z)

@test Graphs.nv(g) == length(vorts_coords)
@test Graphs.ne(g) > 0
@test length(vorts_coords) > 0
@test (length(vort_lines), length(vort_loops), length(vort_rings)) == (16, 0, 1)
@test !isempty(edge_list_periodic)

connected_vorts = VortexDistributions.Detection3D.TimCop.connect_vortex_ends(
    vorts_coords,
    vort_lines,
    vort_loops,
    vort_rings,
    X,
)

@test length(connected_vorts) == 3

g2, lines2, loops2, rings2, coords2, periodic2 = VortexDistributions.full_algorithm(psi, x, y, z)
@test Graphs.nv(g2) == Graphs.nv(g)
@test Graphs.ne(g2) == Graphs.ne(g)
@test length(coords2) == length(vorts_coords)
@test length(lines2) == length(vort_lines)
@test length(loops2) == length(vort_loops)
@test length(rings2) == length(vort_rings)
@test length(periodic2) == length(edge_list_periodic)
