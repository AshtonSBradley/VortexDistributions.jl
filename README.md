# VortexDistributions.jl  

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://AshtonSBradley.github.io/VortexDistributions.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://AshtonSBradley.github.io/VortexDistributions.jl/dev) -->
[![Build Status](https://github.com/AshtonSBradley/VortexDistributions.jl/workflows/CI/badge.svg)](https://github.com/AshtonSBradley/VortexDistributions.jl/actions)
[![Code Quality](https://github.com/AshtonSBradley/VortexDistributions.jl/actions/workflows/Quality.yml/badge.svg)](https://github.com/AshtonSBradley/VortexDistributions.jl/actions/workflows/Quality.yml)
[![TagBot](https://github.com/AshtonSBradley/VortexDistributions.jl/actions/workflows/TagBot.yml/badge.svg)](https://github.com/AshtonSBradley/VortexDistributions.jl/actions/workflows/TagBot.yml)
[![Coverage](https://codecov.io/gh/AshtonSBradley/VortexDistributions.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/AshtonSBradley/VortexDistributions.jl)

Tools for creating and detecting quantum vortices in Bose-Einstein condensates.

## 2D API

The main 2D entry points are:

- `Field`, `Torus`, `Sphere`
- `PointVortex`, `ScalarVortex`
- `findvortices`, `findvortices_grid`, `findvortices_jumps`
- `phase_jumps`, `unwrap`, `zoom_interp`, `zoom_grid`
- `remove_vortices_edge`, `keep_vortices`, `circ_mask`
- `vortex_array` and the current 2D vortex creation helpers

# Installation
```julia
]add VortexDistributions
```

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
## Recursive cluster algorithm

![](/examples/rca_random_example.png)

The recursive cluster algorithm decomposes 2D vortex distributions into
- free vortices
- vortex dipoles
- vortex clusters

using the algorithm developed in 

___Onsager-Kraichnan Condensation in Decaying Two-Dimensional Quantum Turbulence___,\
Thomas P. Billam, Matthew T. Reeves, Brian P. Anderson, and Ashton S. Bradley, \
[Physical Review Letters 112, 145301 (2014)](http://dx.doi.org/10.1103/PhysRevLett.112.145301)

### Usage
```
result = RCA.recursive_cluster_algorithm(vortices, Lx, Ly)
```

For more detail see the usage example script

`/examples/rca_random_example.jl`

## Experimental 3D API

![](/examples/quench_64_subfigures_taller.png)

The 3D detection backend lives under:

```julia
VortexDistributions.Detection3D.TimCop
```

The experimental root entry points are:

```julia
result = detect_vortices_3d(psi, x, y, z)
g, lines, loops, rings, coords, periodic_edges = full_algorithm(psi, x, y, z) # compatibility wrapper
```

The typed `Detection3DResult` returned by `detect_vortices_3d` is the preferred phase 3 API. The older code path is still preserved under `VortexDistributions.Detection3D.Legacy` for comparison only.

The package test suite exercises the vendored backend against the reference fixture at `test/3d/box_vorts.jld2`.

The plotting workflow is kept as runnable example scripts so `GLMakie` and
`GraphMakie` stay out of the main package dependency graph.

### Quick 3D Demo

The quickest end-to-end demonstration uses the bundled `box_vorts.jld2` fixture
and writes a PNG plot of the detected 3D vortex structures. The example
intentionally requires `GLMakie` and `GraphMakie`, but those packages are not
declared as dependencies of `VortexDistributions.jl` itself.

Install the example-only visualization packages in your active environment:

```bash
julia --project=. -e 'using Pkg; Pkg.add(["GLMakie", "GraphMakie"])'
```

Then run the demo:

```bash
julia --project=. examples/timcop_detection_demo.jl
```

Optional arguments:

- first positional argument: output PNG path
- second positional argument: interpolation factor `n_itp`

This example is not part of package testing; it is a runnable usage
demonstration with plotted results, following the original repository’s
`GLMakie` + `GraphMakie` style.

### Full 64^3 Workflow

For the full `quench_64_3d_algo.jl` workflow, install the simulation dependency,
generate the validation dataset, then run either the summary script or the full
plotting script.

1. Install the optional simulation dependency:

```bash
julia --project=. -e 'using Pkg; Pkg.add("FourierGPE")'
```

2. Generate the validation dataset:

```bash
julia --project=. examples/generate_3d_validation_data.jl
```

This writes:

- `examples/3d_validation_data/sim_64.jld2`
- `examples/3d_validation_data/sol_64.jld2`

3. Run the integrated detector and save a summary bundle:

```bash
julia --project=. examples/timcop_deep_test.jl \
  examples/3d_validation_data/sim_64.jld2 \
  examples/3d_validation_data/sol_64.jld2
```

4. Run the full plotting workflow:

```bash
julia --project=. examples/quench_64_3d_algo.jl \
  examples/3d_validation_data/sim_64.jld2 \
  examples/3d_validation_data/sol_64.jld2
```

The plotting script:

- show the iso-surface of the chosen state
- run the 3D detector with `n_itp`
- overlay the detection graph with `GraphMakie`
- apply repeated CCMA smoothing to lines, loops, and rings
- render the final smoothed filaments with `GLMakie`

Optional arguments:

- third argument: simulation time index, default `20`
- fourth argument: interpolation factor `n_itp`, default `4`

The summary script:

- loads the generated 64^3 simulation data
- runs `detect_vortices_3d(psi, x, y, z; n_itp=...)`
- computes connected vortex structures with `Detection3D.connect_vortex_ends`
- prints a structural summary
- saves a JLD2 bundle of the result next to the supplied solution file

This path is intended for exploratory validation and stress testing that is too heavy or dataset-specific for the package test suite.


# External links
___Signatures of Coherent Vortex Structures in a Disordered 2D Quantum Fluid___,\
Matthew T. Reeves, Thomas P. Billam, Brian P. Anderson, and Ashton S. Bradley, \
[Physical Review A 89, 053631 (2014)](http://journals.aps.org/pra/abstract/10.1103/PhysRevA.89.053631)

___Onsager-Kraichnan Condensation in Decaying Two-Dimensional Quantum Turbulence___,\
Thomas P. Billam, Matthew T. Reeves, Brian P. Anderson, and Ashton S. Bradley, \
[Physical Review Letters 112, 145301 (2014)](http://dx.doi.org/10.1103/PhysRevLett.112.145301)
