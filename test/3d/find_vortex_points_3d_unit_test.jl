using JLD2
using Test
using VortexDistributions

@load joinpath(@__DIR__, "box_vorts.jld2") psi_tubes1 X

psi = psi_tubes1
x, y, z = X
tim = VortexDistributions.Detection3D.TimCop

result = VortexDistributions.detect_vortices_3d(psi, x, y, z)
x2, y2, z2, psi2 = tim._interpolate_grid(psi, x, y, z, 2)
interp_result = VortexDistributions.detect_vortices_3d(psi, x, y, z; n_itp = 2)

@test size(psi2) == (128, 128, 128)
@test (length(x2), length(y2), length(z2)) == (128, 128, 128)
@test length(interp_result.coords) > length(result.coords)
@test !isempty(interp_result.periodic_edges)
@test length(interp_result.vortex_rings) == 1
