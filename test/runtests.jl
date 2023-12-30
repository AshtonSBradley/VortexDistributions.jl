using VortexDistributions, Test, SafeTestsets

@safetestset "No Vortex Fallback" begin include("no_vortex_test.jl") end

@safetestset "Point Vortex" begin include("point_vortex_test.jl") end

@safetestset "Single Vortex Accuracy" begin include("single_vort_acc_test.jl") end

@safetestset "Multiple Vortices" begin include("multivort_test.jl") end

@safetestset "Periodic Dipole" begin include("periodic_dipole_test.jl") end 

@safetestset "Creation functions" begin include("creations_test.jl") end 

# @safetestset "Multiple Vortices 2" begin include("multivort_test2.jl") end
@safetestset "Unwrap" begin include("test_unwrap.jl") end

@safetestset "Unwrap2D" begin include("test_unwrap2D.jl") end

@safetestset "Phase Jumps" begin include("phase_jumps.jl") end

@safetestset "Scalar Vortex Arguments" begin include("vortex_types.jl") end

@safetestset "Single Vortex Sphere Field" begin include("single_vort_sphere.jl") end

@safetestset "Keep Vortices" begin include("keep_vortices_test.jl") end

@safetestset "Vortex 3D example" begin include("3d/box_vorts_example_test.jl") end

@safetestset "find_vortex_points_3d unit tests" begin include("3d/find_vortex_points_3d_unit_test.jl") end

@safetestset "connect_vortex_points_3d unit tests" begin include("3d/connect_vortex_points_3d_unit_tests.jl") end

@safetestset "sort_classified_vorts_3d unit tests" begin include("3d/sort_classified_vorts_3d_unit_tests.jl") end