using JLD2
using Test
using VortexDistributions

@load joinpath(@__DIR__, "box_vorts.jld2") psi_tubes1 X

psi = psi_tubes1
x, y, z = X

result = VortexDistributions.detect_vortices_3d(psi, x, y, z)
connected_from_result = VortexDistributions.Detection3D.connect_vortex_ends(result, X)
connected_from_parts = VortexDistributions.Detection3D.connect_vortex_ends(
    result.coords,
    result.vortex_lines,
    result.vortex_loops,
    result.vortex_rings,
    X,
)

@test connected_from_result == connected_from_parts
@test length(connected_from_result) == 3
@test length(last(connected_from_result)) == 1

singleton = VortexDistributions.Detection3D.connect_vortex_ends(result.coords, [result.vortex_lines[1]], [], [], X)
@test singleton == [[result.vortex_lines[1]]]
