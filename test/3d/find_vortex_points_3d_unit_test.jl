using VortexDistributions, JLD2, Test

@load joinpath(@__DIR__, "box_vorts.jld2") psi_tubes1 X

psi = psi_tubes1


@test_throws AssertionError find_vortex_points_3d(psi[:, :, 1:10], X, 1)

@test_throws AssertionError find_vortex_points_3d(psi[1:10, :, :], X, 1)

@test_throws AssertionError find_vortex_points_3d(psi, X, 17)

@test_throws MethodError find_vortex_points_3d(psi[1, :, :], X)

@test_throws MethodError find_vortex_points_3d(psi, (X[1], X[2]))

@test_throws MethodError find_vortex_points_3d(psi, X, 1.1)

@test typeof(find_vortex_points_3d(psi, X, 1)) ==  Vector{Vector{Float64}}