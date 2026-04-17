# Recursive Cluster Algorithm Status

## Current state

The recursive cluster algorithm (RCA) is now part of the package proper and
lives under [`src/RecursiveClusterAlgorithm.jl`](../src/RecursiveClusterAlgorithm.jl)
with helper files in [`src/RecursiveClusterAlgorithm/`](../src/RecursiveClusterAlgorithm).

The supported public entry points are:

- `VortexDistributions.RecursiveClusterAlgorithm`
- `VortexDistributions.RCAResult`
- `VortexDistributions.recursive_cluster_algorithm`

The package reuses shared vortex primitives from [`src/Core/`](../src/Core),
including `PointVortex`, `Dipole`, `Cluster`, and the internal periodic-distance
and spanning-tree helpers.

## Test coverage

RCA behavior is covered by the main package test suite rather than a separate
prototype test tree.

Relevant coverage lives in:

- [`test/recursive_cluster_algorithm_test.jl`](../test/recursive_cluster_algorithm_test.jl)
- [`test/runtests.jl`](../test/runtests.jl)

The current tests exercise:

- dipole extraction
- seed selection
- positive and negative cluster growth
- the top-level RCA driver
- edge cases such as all-dipole inputs and invalid array lengths

## Cleanup note

The old `_research/` prototype tree has been removed from the repository.
That code was no longer part of the package load path, duplicated shipped RCA
logic, and contained unfinished MATLAB-era experiments and one-off scripts.

## Future work

The remaining RCA improvements are package work, not migration work. The most
useful next steps are:

- add more periodic-boundary regression cases
- add tests for larger mixed configurations and classification stability
- decide whether RCA needs richer result metadata without changing the core
  `Cluster` type
