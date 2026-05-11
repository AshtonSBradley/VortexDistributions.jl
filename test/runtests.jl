using Test
using VortexDistributions

@testset "No Vortex Fallback" begin
    include("no_vortex_test.jl")
end

@testset "Point Vortex" begin
    include("point_vortex_test.jl")
end

@testset "Single Vortex Accuracy" begin
    include("single_vort_acc_test.jl")
end

@testset "Multiple Vortices" begin
    include("multivort_test.jl")
end

@testset "Detection Helpers" begin
    include("multivort_test2.jl")
end

@testset "Periodic Dipole" begin
    include("periodic_dipole_test.jl")
end

@testset "Creation functions" begin
    include("creations_test.jl")
end

@testset "Unwrap" begin
    include("test_unwrap.jl")
end

@testset "Unwrap2D" begin
    include("test_unwrap2D.jl")
end

@testset "Interpolation" begin
    include("test_interp.jl")
end

@testset "Phase Jumps" begin
    include("phase_jumps.jl")
end

@testset "Scalar Vortex Arguments" begin
    include("vortex_types.jl")
end

@testset "Single Vortex Sphere Field" begin
    include("single_vort_sphere.jl")
end

@testset "Keep Vortices" begin
    include("keep_vortices_test.jl")
end

@testset "Ripley Clustering" begin
    include("ripley_test.jl")
end

@testset "API Regression" begin
    include("api_regression_test.jl")
end

@testset "Core Assets" begin
    include("new_jld_files.jl")
end

@testset "Recursive Cluster Algorithm" begin
    include("recursive_cluster_algorithm_test.jl")
end

@testset "3D TimCop Helpers" begin
    include("3d/box_vorts_example_test.jl")
end

@testset "3D Connection Helpers" begin
    include("3d/connect_vortex_points_3d_unit_tests.jl")
end

@testset "3D Interpolation" begin
    include("3d/find_vortex_points_3d_unit_test.jl")
end

@testset "3D TimCop Ring Handling" begin
    include("3d/sort_classified_vorts_3d_unit_tests.jl")
end

@testset "TimCop 3D Backend" begin
    include("3d/timcop_backend_test.jl")
end
