# VortexDistributions.jl <img align="right" src="/examples/vortfluid.gif" width="100" height="100">

[![Build Status](https://travis-ci.org/AshtonSBradley/VortexDistributions.jl.svg?branch=master)](https://travis-ci.org/AshtonSBradley/VortexDistributions.jl)  [![Coverage Status](https://coveralls.io/repos/AshtonSBradley/VortexDistributions.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/AshtonSBradley/VortexDistributions.jl?branch=master)  [![codecov.io](http://codecov.io/github/AshtonSBradley/VortexDistributions.jl/coverage.svg?branch=master)](http://codecov.io/github/AshtonSBradley/VortexDistributions.jl?branch=master)

Tools for working with distributions of two-dimensional quantum vortices in Bose-Einstein condendates.

- [x] Fast, accurate vortex detection.
  - Uses a highly optimized version of the plaquette method (phase intergral around each 4-point plaquette), with recursive interpolation to acheive a good balance between speed and accuracy.
  - At present only tests for charge +/-1 in 2D

- [x] Vortex creation
  - Solves the 2D GPE problem for charge n on the infinite domain
  - Interpolates vortex solution to density and phase imprint on arbitrary 2D domains
- [ ] Recursive cluster algorithm
- [ ] Vortex correlation functions

# Detection Example
```julia
using VortexDistributions, Plots
gr(transpose=true,xlabel="x",ylabel="y",legend=false)

# Our example system has grid point spacing
# of two points per healing length
# (the following is in units of healing length)
Lx=200;Nx=400;
Ly=200;Ny=400
x = linspace(-Lx/2,Lx/2,Nx);dx=diff(x)[1]
y = linspace(-Ly/2,Ly/2,Ny);dy=diff(y)[1]

#make charge-1 vortex near the point (x,y)=(10,10)
facx,facy = rand(2)
testvort = [10+dx*facx 10+dy*facy 1.0]

#construct vortex wavefunction, including both density and phase
ψ = one.(x.*y') |> complex
makevortex!(ψ,testvort,x,y);
```

In this example the vortex is created at
```julia
julia> testvort
Out[25]:
1×3 Array{Float64,2}:
 10.3234  10.314  1.0
 ```
 We can find all the vortices using (removing edge vortices by default):
 ```julia
 nt,np,nn,vortices = findvortices(ψ,x,y)
 ```
 For our single vortex example, the vortex is detected with `(x,y,charge)`
 ```julia
 julia> vortices
Out[26]:
1×3 Array{Float64,2}:
 10.3274  10.3178  1.0
 ```

Plotting the results, we have the phase at successive zoom levels with vortex location, `+`, and detected location, `o` (see examples):

![](/examples/phase.png)

and density at successive zoom levels with vortex location and detected location:

![](/examples/density.png)

 The benchmark gives (2015 MacBook Air 1.6GHz i5)
 ```julia
 using BenchmarkTools
 julia> @btime nt,np,nn,vortices = findvortices(ψ,x,y)
  6.165 ms (550 allocations: 3.96 MiB)
 ```

#### Acknowledgements
Matthew Reeves, Thomas Billam, Michael Cawte

# External links
___Signatures of Coherent Vortex Structures in a Disordered 2D Quantum Fluid___,\
Matthew T. Reeves, Thomas P. Billam, Brian P. Anderson, and Ashton S. Bradley, \
[Physical Review A 89, 053631 (2014)](http://journals.aps.org/pra/abstract/10.1103/PhysRevA.89.053631)

___Onsager-Kraichnan Condensation in Decaying Two-Dimensional Quantum Turbulence___,\
Thomas P. Billam, Matthew T. Reeves, Brian P. Anderson, and Ashton S. Bradley, \
[Physical Review Letters 112, 145301 (2014)](http://dx.doi.org/10.1103/PhysRevLett.112.145301)
