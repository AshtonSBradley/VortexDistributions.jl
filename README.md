# VortexDistributions.jl  

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://AshtonSBradley.github.io/VortexDistributions.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://AshtonSBradley.github.io/VortexDistributions.jl/dev) -->
[![Build Status](https://github.com/AshtonSBradley/VortexDistributions.jl/workflows/CI/badge.svg)](https://github.com/AshtonSBradley/VortexDistributions.jl/actions)
[![Coverage](https://codecov.io/gh/AshtonSBradley/VortexDistributions.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/AshtonSBradley/VortexDistributions.jl)

Tools for creating and detecting quantum vortices in Bose-Einstein condensates.

## Stable 2D API

Phase 1 refactoring keeps the established 2D user-facing API at the package root:

- `Field`, `Torus`, `Sphere`
- `PointVortex`, `ScalarVortex`
- `findvortices`, `findvortices_grid`, `findvortices_jumps`
- `phase_jumps`, `unwrap`, `zoom_interp`, `zoom_grid`
- `remove_vortices_edge`, `keep_vortices`, `circ_mask`
- `vortex_array` and the current 2D vortex creation helpers

Internally, the package is now split into:

- `VortexDistributions.Core`
- `VortexDistributions.Detection2D`
- `VortexDistributions.Creation2D`
- `VortexDistributions.Analysis2D`
- `VortexDistributions.Detection3D.Legacy`

Low-level helpers are no longer exported from the root namespace.

## Experimental 3D API

Current 3D code is treated as experimental / legacy and is intentionally isolated under:

```julia
VortexDistributions.Detection3D.Legacy
```

The legacy entry point `VortexDistributions.full_algorithm` is still available as a compatibility shim, but new work should target the experimental namespace until the phase 2 backend lands.

# Installation
```julia
]add VortexDistributions
```

See also:

- [Phase 1 migration note](docs/phase1-migration.md)
- [Phase 2 integration map](docs/phase2-integration-map.md)

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
