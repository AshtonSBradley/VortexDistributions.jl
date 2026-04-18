using Graphs
using JLD2
using Test
using VortexDistributions

@load joinpath(@__DIR__, "box_vorts.jld2") psi_tubes1 X

psi = psi_tubes1
x, y, z = X

legacy_backend = VortexDistributions.Detection3D.LegacyPlaneBackend3D()
legacy_result = VortexDistributions.detect_vortices_3d(psi, x, y, z; backend = legacy_backend)
legacy_full = VortexDistributions.full_algorithm(psi, x, y, z; backend = legacy_backend)

@test legacy_result isa VortexDistributions.Detection3D.Detection3DResult
@test legacy_result.periodic_edges == Tuple{Int, Int}[]
@test (length(legacy_result.vortex_lines), length(legacy_result.vortex_loops), length(legacy_result.vortex_rings)) == (16, 0, 1)
@test length(legacy_result.coords) == 958

g, lines, loops, rings, coords, periodic_edges = legacy_full
@test Graphs.nv(g) == 958
@test Graphs.ne(g) == 942
@test length(lines) == 16
@test isempty(loops)
@test length(rings) == 1
@test length(coords) == 958
@test isempty(periodic_edges)
