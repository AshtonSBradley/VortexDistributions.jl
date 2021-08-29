using VortexDistributions

## Basic unit testing for coverage 
@test scalar_ansatz(1) == scalar_ansatz(-1)
f = Ansatz()
x = LinRange(0, 10, 100)
y = f.(x)
@test length(y) == length(x)
@test f.(3,4) == f(5)

y, ψ, res = VortexDistributions.gpecore_exact(1,2,100)
@test length(y) == 100