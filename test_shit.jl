using VortexDistributions, Parameters, LinearAlgebra

foo = Ansatz()
@unpack f, ξ, Λ = foo
x = LinRange(0, 10, 100)
y = f.(x)

foo.(3,4)
foo(5)

using VortexDistributions
vorts = found_near(30, 10)

