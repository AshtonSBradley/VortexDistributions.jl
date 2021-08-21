using VortexDistributions, Test

# @testset "Creation Functions" begin 
@test scalar_ansatz(0) == 0
@test scalar_ansatz(1) == scalar_ansatz(-1)
@test typeof(Ansatz()) <: Ansatz 
# end 
