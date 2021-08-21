using VortexDistributions, Test

# @testset "Creation Functions" begin 
@test scalar_ansatz(0) == 0
@test scalar_ansatz(1) == scalar_ansatz(-1)
f = Ansatz()
x = LinRange(0, 10, 100)
y = f.(x)
@test length(y) == length(x)
@test f.(3,4) == f(5)
