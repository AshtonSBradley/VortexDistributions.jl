# Phase 1 Migration Note

`VortexDistributions.jl` now has an internal module layout aimed at keeping the stable 2D API intact while making room for a future 3D backend integration.

## Stable 2D API

These root-level APIs are still the supported package entry points:

- `Field`, `Torus`, `Sphere`
- `PointVortex`, `ScalarVortex`
- `findvortices`, `findvortices_grid`, `findvortices_jumps`
- `phase_jumps`, `unwrap`, `zoom_interp`, `zoom_grid`
- `remove_vortices_edge`, `keep_vortices`, `circ_mask`
- `vortex_array` and the current 2D vortex creation helpers

## User-visible changes

- Internal code is now split across `Core`, `Detection2D`, `Creation2D`, `Analysis2D`, and `Detection3D.Legacy`.
- Low-level helper functions such as `Δ` and `find_where` are no longer exported from the package root.
- The current 3D path is explicitly experimental and lives under `VortexDistributions.Detection3D.Legacy`.
- `VortexDistributions.full_algorithm` remains available as a compatibility shim for the legacy 3D entry point.

## Dependency hygiene

- Heavy or stale imports that were not on the stable 2D path have been removed from the main package load path.
- The new backend abstraction points isolate future detection-engine swaps from the public 2D API.
