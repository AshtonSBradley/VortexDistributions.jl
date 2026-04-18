using Graphs
using JLD2
using Test
using VortexDistributions

@load joinpath(@__DIR__, "box_vorts.jld2") psi_tubes1 X

psi = psi_tubes1
x, y, z = X
tim = VortexDistributions.Detection3D.TimCop

Δϕx, Δϕy, Δϕz = tim.findvortices_planes_threaded(psi)
edge_list, edge_list_periodic, vorts_x, vorts_y, vorts_z = tim.vortex_edge_list_threaded(Δϕx, Δϕy, Δϕz)
coords = tim.vortex_coords(vorts_x, vorts_y, vorts_z, x, y, z)
g = SimpleGraph(Edge.(edge_list))
vort_lines, vort_loops, vort_rings = tim.link_graph_vorts(g)

@test (length(vorts_x), length(vorts_y), length(vorts_z)) == (364, 308, 286)
@test length(edge_list) == 1884
@test !isempty(edge_list_periodic)
@test length(coords) == 958
@test Graphs.nv(g) == 958
@test Graphs.ne(g) == 942
@test (length(vort_lines), length(vort_loops), length(vort_rings)) == (16, 0, 1)
