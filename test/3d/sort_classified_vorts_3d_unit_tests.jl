using JLD2
using Test
using VortexDistributions

@load joinpath(@__DIR__, "box_vorts.jld2") psi_tubes1 X

psi = psi_tubes1
x, y, z = X
tim = VortexDistributions.Detection3D.TimCop

result = VortexDistributions.detect_vortices_3d(psi, x, y, z)
connected = VortexDistributions.Detection3D.connect_vortex_ends(result, X)

_, lines_no_repeat, loops_no_repeat, rings_no_repeat, coords_no_repeat, periodic_no_repeat =
    tim.full_algorithm(psi, x, y, z; repeat_end_of_ring = false)

@test length(connected) == 3
@test sort(length.(connected)) == [1, 4, 12]

@test length(result.vortex_rings) == 1
@test length(rings_no_repeat) == 1
@test length(result.vortex_rings[1]) == length(rings_no_repeat[1]) + 1
@test first(result.vortex_rings[1]) == last(result.vortex_rings[1])
@test first(rings_no_repeat[1]) != last(rings_no_repeat[1])

@test length(lines_no_repeat) == length(result.vortex_lines)
@test isempty(loops_no_repeat)
@test length(coords_no_repeat) == length(result.coords)
@test length(periodic_no_repeat) == length(result.periodic_edges)
