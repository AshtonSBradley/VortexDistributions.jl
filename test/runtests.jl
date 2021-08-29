using VortexDistributions, Test, SafeTestsets

@safetestset "Point Vortex" begin include("point_vortex_test.jl") end
@safetestset "Single Vortex Accuracy" begin include("single_vort_acc_test.jl") end
@safetestset "Multiple Vortices" begin include("multivort_test.jl") end
@safetestset "Periodic Dipole" begin include("periodic_dipole_test.jl") end 
@safetestset "Creation functions" begin include("creations_test.jl") end 
@safetestset "Multiple Vortices 2" begin include("multivort_test2.jl") end
@safetestset "Unwrap" begin include("test_unwrap.jl") end
@safetestset "Unwrap2D" begin include("test_unwrap2D.jl") end
@safetestset "Phase Jumps" begin include("phase_jumps.jl") end