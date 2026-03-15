# Phase 2 Integration Map

This note describes where `VortexDetection.jl` components should land when phase 2 starts.

## 2D integration boundary

- `src/Detection2D/backends.jl`
- `src/Detection2D/detection.jl`

`AbstractDetectionBackend2D` and `PlaquetteBackend2D` are the current seam. A future `VortexDetection2DBackend` should implement the candidate-finding methods used by `findvortices_grid` and `findvortices_jumps` without changing the root 2D API.

## 3D integration boundary

- `src/Detection3D.jl`
- `src/Detection3D/legacy.jl`

`Detection3D.Legacy` is the quarantine zone for the old threaded plane-linking approach. Phase 2 should add a new backend alongside it first, keep the legacy implementation available during transition, and only promote a new root-facing 3D API after validation.

## Suggested landing zones

- New backend interfaces and shared traits: `src/Detection2D/backends.jl`, `src/Detection3D.jl`
- Shared geometry or classification helpers reused by 2D and 3D backends: `src/Core/`
- Stable public adapters that preserve current user entry points: `src/VortexDistributions.jl`

## TODO markers to follow

- `src/Detection2D/backends.jl`
- `src/Detection3D.jl`

Those files contain the phase 2 TODO anchors added during phase 1.
