# Recursive Cluster Algorithm Merge Plan

## Goal

Merge the usable recursive cluster algorithm (RCA) code from [`_research/`](../_research) into the package without rewriting the established `src` core. The only planned `src` changes are the minimum needed to expose RCA as a supported surface.

## Constraints

- Do not rewrite the existing `src/Core`, `src/Creation2D`, `src/Analysis2D`, `src/Detection2D`, or `src/Detection3D` implementations.
- Reuse the primitives that already exist in `src/Core`.
- Remove duplicate RCA definitions that are still living in `_research`.
- Keep unfinished MATLAB-style ports and one-off scripts quarantined until they are translated and tested.

## What Already Exists In `src`

The package already absorbed several RCA building blocks:

- [`src/Core/types.jl`](../src/Core/types.jl): `Cluster`, `Dipole`, `VortexGroup`
- [`src/Core/internal.jl`](../src/Core/internal.jl): `distances`, `periodic_distances!`, `sparse_distances`, `spanning_tree`

That means `_research/RecursiveClusterAlgorithm` should not continue owning those definitions.

## Duplication To Remove

### Direct duplicates

- [`_research/RecursiveClusterAlgorithm/RecursiveClusterAlgorithm.jl`](../_research/RecursiveClusterAlgorithm/RecursiveClusterAlgorithm.jl)
  - redefines `Cluster`
  - redefines the `Cluster(vort)` constructor
  - re-exports graph/distance helpers that already live in `src/Core`
- [`_research/RecursiveClusterAlgorithm/utils.jl`](../_research/RecursiveClusterAlgorithm/utils.jl)
  - duplicates `distances`
  - duplicates `periodic_distances!`
  - duplicates `make_periodic!`
  - duplicates `sparse_distances`
  - duplicates `spanning_tree`
- Root-level prototype files in [`_research/`](../_research)
  - `get_dipoles_rca.jl`
  - `seed_clusters_rca.jl`
  - `grow_plus_clusters_rca.jl`
  - `grow_minus_clusters_rca.jl`
  - these are older parallel copies of the same RCA work and should not survive the merge

### Code not ready for promotion

These files are still incomplete or still effectively MATLAB syntax:

- [`_research/RecursiveClusterAlgorithm/RCA.jl`](../_research/RecursiveClusterAlgorithm/RCA.jl)
- [`_research/RecursiveClusterAlgorithm/get_spanning_trees.jl`](../_research/RecursiveClusterAlgorithm/get_spanning_trees.jl)
- [`_research/RecursiveClusterAlgorithm/get_positive_spanning_trees.jl`](../_research/RecursiveClusterAlgorithm/get_positive_spanning_trees.jl)
- [`_research/RecursiveClusterAlgorithm/get_negative_spanning_trees.jl`](../_research/RecursiveClusterAlgorithm/get_negative_spanning_trees.jl)
- [`_research/RecursiveClusterAlgorithm/find_smallest_spanning_trees.jl`](../_research/RecursiveClusterAlgorithm/find_smallest_spanning_trees.jl)
- [`_research/RecursiveClusterAlgorithm/sort_RCA_structs.jl`](../_research/RecursiveClusterAlgorithm/sort_RCA_structs.jl)

## Recommended Landing Shape

### 1. Keep the algorithm implementation sourced from `_research`

Do not inline RCA logic into `src/Core`. Instead, treat `_research/RecursiveClusterAlgorithm` as the algorithm source of truth during this merge.

### 2. Add a thin wrapper module in `src`

Add one new wrapper file under `src`, for example `src/RecursiveClusterAlgorithm.jl`, that:

- defines a package-owned submodule, for example `module RecursiveClusterAlgorithm`
- imports the existing package primitives from `..Core`
- includes the translated RCA implementation files from `_research/RecursiveClusterAlgorithm`
- exposes only the supported RCA entry points

This keeps `src` changes limited to exposure, not core rewrites.

### 3. Extend the package root with one include and one export block

In [`src/VortexDistributions.jl`](../src/VortexDistributions.jl):

- `include("RecursiveClusterAlgorithm.jl")`
- bind the submodule as a stable namespace, for example `const RecursiveClusterAlgorithm = RecursiveClusterAlgorithm`
- optionally re-export only the top-level RCA entry point if a flat API is desired

Prefer keeping RCA primarily namespaced under `VortexDistributions.RecursiveClusterAlgorithm` to avoid crowding the root API.

## Merge Phases

### Phase 1: Stabilize the minimal RCA surface

Promote only the pieces that are already close to Julia-ready:

- [`_research/RecursiveClusterAlgorithm/get_dipoles.jl`](../_research/RecursiveClusterAlgorithm/get_dipoles.jl)
- [`_research/RecursiveClusterAlgorithm/seed_clusters.jl`](../_research/RecursiveClusterAlgorithm/seed_clusters.jl)
- [`_research/RecursiveClusterAlgorithm/grow_plus_clusters.jl`](../_research/RecursiveClusterAlgorithm/grow_plus_clusters.jl)
- [`_research/RecursiveClusterAlgorithm/grow_minus_clusters.jl`](../_research/RecursiveClusterAlgorithm/grow_minus_clusters.jl)

Before exposing them:

- rewrite them to call `Core.distances`, `Core.periodic_distances!`, and `Core.spanning_tree` instead of local duplicates
- normalize return shapes and types
- replace MATLAB-style scalar logic and indexing assumptions
- remove plotting dependencies from the core algorithm path

Recommended public surface for phase 1:

- `RecursiveClusterAlgorithm.get_dipoles`
- `RecursiveClusterAlgorithm.seed_clusters`
- `RecursiveClusterAlgorithm.grow_plus_clusters`
- `RecursiveClusterAlgorithm.grow_minus_clusters`

### Phase 2: Add one package-level RCA driver

After phase 1 helpers are tested, add a single orchestrating entry point, for example:

- `RecursiveClusterAlgorithm.recursive_cluster_algorithm`

That driver should live in `_research` implementation files, but be loaded through the `src` wrapper. It should:

- accept `Vector{PointVortex}` or separate coordinate/charge arrays
- recursively identify dipoles
- seed same-sign clusters
- grow positive and negative clusters
- return package-native `Dipole`, `Cluster`, and leftover free vortices

Do not port the current [`RCA.jl`](../_research/RecursiveClusterAlgorithm/RCA.jl) file verbatim. Use it only as algorithm reference, because it still contains MATLAB control flow, `eval`, `disp`, and struct-building patterns that do not fit the package.

### Phase 3: Periodic spanning tree support

Only after the RCA driver is working:

- translate periodic spanning-tree minimization from the incomplete `_research` files
- represent cluster metadata with Julia structs or named tuples, not dynamic `cluster1`, `cluster2`, ... fields
- fold the result back into package-native `Cluster` values instead of parallel RCA-only cluster containers

This phase should also decide whether `Cluster` needs extra metadata. If yes, add it as an RCA result wrapper rather than changing the existing core `Cluster` type.

## Concrete File Actions

### Keep and translate

- `_research/RecursiveClusterAlgorithm/get_dipoles.jl`
- `_research/RecursiveClusterAlgorithm/seed_clusters.jl`
- `_research/RecursiveClusterAlgorithm/grow_plus_clusters.jl`
- `_research/RecursiveClusterAlgorithm/grow_minus_clusters.jl`

### Keep as reference only for now

- `_research/RecursiveClusterAlgorithm/RCA.jl`
- `_research/RecursiveClusterAlgorithm/get_spanning_trees.jl`
- `_research/RecursiveClusterAlgorithm/get_positive_spanning_trees.jl`
- `_research/RecursiveClusterAlgorithm/get_negative_spanning_trees.jl`
- `_research/RecursiveClusterAlgorithm/find_smallest_spanning_trees.jl`
- `_research/RecursiveClusterAlgorithm/sort_RCA_structs.jl`

### Delete after wrapper is in place

- `_research/RecursiveClusterAlgorithm/utils.jl`
- `_research/get_dipoles_rca.jl`
- `_research/seed_clusters_rca.jl`
- `_research/grow_plus_clusters_rca.jl`
- `_research/grow_minus_clusters_rca.jl`

### Leave out of the package load path

- plotting helpers under `_research/RecursiveClusterAlgorithm/plot_*.jl`
- ad hoc scripts in `_research/*.jl`

Those can stay as experiments/examples until there is a clear package API for visualization.

## Testing Plan

Add RCA-focused tests to the main package test suite instead of relying on `_research/RecursiveClusterAlgorithm/test`.

Recommended new test groups:

- `test/recursive_cluster/get_dipoles_test.jl`
  - empty input
  - one mutual opposite-sign pair
  - periodic dipole across boundaries
- `test/recursive_cluster/seed_clusters_test.jl`
  - mutual same-sign nearest neighbors
  - mixed-sign configurations that should not seed
- `test/recursive_cluster/grow_clusters_test.jl`
  - positive cluster growth
  - negative cluster growth
  - no-growth edge case
- `test/recursive_cluster/api_test.jl`
  - package namespace exposure
  - result types are `Dipole`, `Cluster`, and `PointVortex`

Integrate them from [`test/runtests.jl`](../test/runtests.jl) only after the wrapper module is wired in.

## Recommended Implementation Order

1. Add the `src` wrapper module and root include/export surface.
2. Repoint `_research/RecursiveClusterAlgorithm` code to `src/Core` primitives.
3. Delete duplicate RCA copies in `_research`.
4. Translate and test `get_dipoles` and `seed_clusters`.
5. Translate and test `grow_plus_clusters` and `grow_minus_clusters`.
6. Add the single RCA driver entry point.
7. Tackle periodic spanning-tree refinement as a later pass.

## Decision Summary

The cleanest merge is not to move `_research` wholesale into `src/Core`. The repo is already set up with the right split: `src/Core` owns shared vortex and graph primitives, while RCA should be exposed as a separate higher-level algorithm surface that reuses those primitives. That removes duplication, respects the current `src` boundary, and avoids importing unfinished MATLAB-era code into the stable package path.
