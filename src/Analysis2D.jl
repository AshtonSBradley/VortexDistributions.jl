module Analysis2D

using Interpolations

using ..Core: Field, PointVortex, pos, Δ

include("Analysis2D/phase_tools.jl")

export circ_mask, keep_vortices, phase_jumps, phase_jumps!, unwrap, unwrap!, zoom_grid, zoom_interp

end
