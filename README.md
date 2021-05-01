# VortexDistributions.jl <img align="right" src="/examples/vortfluid.gif" width="100" height="100">

<!-- [![Build Status](https://travis-ci.org/AshtonSBradley/VortexDistributions.jl.svg?branch=master)](https://travis-ci.org/AshtonSBradley/VortexDistributions.jl)  [![Coverage Status](https://coveralls.io/repos/AshtonSBradley/VortexDistributions.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/AshtonSBradley/VortexDistributions.jl?branch=master)  [![codecov.io](http://codecov.io/github/AshtonSBradley/VortexDistributions.jl/coverage.svg?branch=master)](http://codecov.io/github/AshtonSBradley/VortexDistributions.jl?branch=master) -->

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://AshtonSBradley.github.io/VortexDistributions.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://AshtonSBradley.github.io/VortexDistributions.jl/dev)
[![Build Status](https://github.com/AshtonSBradley/VortexDistributions.jl/workflows/CI/badge.svg)](https://github.com/AshtonSBradley/VortexDistributions.jl/actions)
[![Coverage](https://codecov.io/gh/AshtonSBradley/VortexDistributions.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/AshtonSBradley/VortexDistributions.jl)

Tools for creating and detecting quantum vortices in Bose-Einstein condensates.
- [x] Fast, accurate vortex detection.
  - Highly optimized version of the plaquette method (phase integral around each 4-point plaquette), with recursive interpolation to achieve a good balance between speed and accuracy.
  - At present only tests for charge +/-1 in 2D
- [x] Vortex creation
  - Solves the 2D GPE problem for charge n on the infinite domain
  - Interpolates vortex solution to density and phase imprint on arbitrary 2D domains
- [ ] Recursive cluster algorithm
- [ ] Vortex correlation functions

# Detection Example
```julia
using VortexDistributions, Plots
gr(xlabel="x",ylabel="y",legend=false)

# make a simple 2D test field
Nx = 400; Ny = Nx
Lx = 200; Ly = Lx
x = LinRange(-Lx / 2, Ly / 2, Nx); y = x
psi0 = one.(x*y') |> complex

# doubly periodic boundary conditions
psi = Torus(psi0,x,y)

# make a point vortex
pv = PointVortex(30.0,70.3,-1)

# make a scalar GPE vortex with exact core
spv = ScalarVortex(pv)
vortex!(psi,spv)

# make some more random vortices
vort = rand_vortex(10,psi)
vortex!(psi,vort)
```

We can recover the raw point vortex data from `PointVortex()` with
```julia
vortex_array(pv)
 ```
 or from a `ScalarVortex()` with
 ```julia
vortex_array(spv.vort)
  ```
 We can find all the vortices, removing edge vortices by default:
 ```julia
vfound = findvortices(psi)
 ```

For a single vortex example, we show have the phase at successive zoom levels with vortex location, `+`, and detected location, `o` (see examples):

![](/examples/phase.png)

and density at successive zoom levels with vortex location and detected location:

![](/examples/density.png)

 The benchmark gives (2018 MacBook Pro 2.33GHz Intel i5)
 ```julia
 using BenchmarkTools
 julia> @btime vort = findvortices(psi)
   4.037 ms (585 allocations: 3.84 MiB)
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
