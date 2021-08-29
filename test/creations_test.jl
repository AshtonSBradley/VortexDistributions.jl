using VortexDistributions

## Basic unit testing for coverage 
@test scalar_ansatz(1) == scalar_ansatz(-1)
f = Ansatz()
x = LinRange(0, 10, 100)
y = f.(x)
@test length(y) == length(x)
@test f.(3,4) == f(5)

# @test VortexDistributions.H(-1) == 0.0
# @test VortexDistributions.H(0.1) == 1.0

# @test VortexDistributions.shift(1, 1) == 0
# @test VortexDistributions.shift(2, 1) == 1

# @test VortexDistributions.tans(1, 1) == tan((0 - π)*0.5)
# @test VortexDistributions.tanhs(1,1,2) == tanh((0 + 2*π*2)*0.5)

# VortexDistributions.kernel(1,1,2,2,3,3,4) == atan()