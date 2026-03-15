using Test
using VortexDistributions

@test isdefined(VortexDistributions, :Core)
@test isdefined(VortexDistributions, :Detection2D)
@test isdefined(VortexDistributions, :Creation2D)
@test isdefined(VortexDistributions, :Analysis2D)
@test isdefined(VortexDistributions, :Detection3D)
@test isdefined(VortexDistributions, :Detection3DLegacy)
@test isdefined(VortexDistributions, :detect_vortices_3d)

for name in (
    :Field,
    :Torus,
    :Sphere,
    :PointVortex,
    :ScalarVortex,
    :findvortices,
    :findvortices_grid,
    :findvortices_jumps,
    :phase_jumps,
    :unwrap,
    :zoom_interp,
    :zoom_grid,
    :remove_vortices_edge,
    :keep_vortices,
    :circ_mask,
    :vortex_array,
    :scalar_ansatz,
    :vortex!,
    :dipole_phase,
    :periodic_dipole!,
    :rand_charge,
    :rand_pointvortex,
    :rand_scalarvortex,
    :rand_vortex,
    :rand_vortexfield,
)
    @test isdefined(VortexDistributions, name)
end

@test !(:Δ in names(VortexDistributions))
@test !(:find_where in names(VortexDistributions))
@test VortexDistributions.full_algorithm === VortexDistributions.Detection3D.full_algorithm
