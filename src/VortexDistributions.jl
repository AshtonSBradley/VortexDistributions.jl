module VortexDistributions

include("Core.jl")
include("Analysis2D.jl")
include("Creation2D.jl")
include("Detection2D.jl")
include("Detection3D.jl")
include("RecursiveClusterAlgorithm.jl")

using .Core: Ansatz,
    ChannelVortex,
    Cluster,
    CoreShape,
    Dipole,
    Exact,
    Field,
    PointVortex,
    ScalarVortex,
    Sphere,
    Torus,
    Vortex,
    VortexGroup,
    charge,
    pos,
    xpos,
    ypos,
    vortex_array
using .Analysis2D: circ_mask, keep_vortices, phase_jumps, phase_jumps!, unwrap, unwrap!, zoom_grid, zoom_interp
using .Creation2D: gpecore_exact,
    periodic_dipole!,
    rand_charge,
    rand_scalarvortex,
    rand_vortex,
    rand_vortexfield,
    scalar_ansatz,
    thetad,
    dipole_phase,
    vortex!
using .Detection2D: findvortices, findvortices_grid, findvortices_jumps, found_near, remove_vortices_edge
using .Detection3D: detect_vortices_3d, full_algorithm
using .RecursiveClusterAlgorithm: RCAResult, recursive_cluster_algorithm

const Detection3DLegacy = Detection3D.Legacy

export Field, Torus, Sphere
export PointVortex, ScalarVortex
export findvortices, findvortices_grid, findvortices_jumps
export phase_jumps, phase_jumps!, unwrap, unwrap!, zoom_interp, zoom_grid
export remove_vortices_edge, keep_vortices, circ_mask
export vortex_array
export scalar_ansatz, vortex!, dipole_phase, periodic_dipole!
export rand_charge, rand_pointvortex, rand_scalarvortex, rand_vortex, rand_vortexfield
export thetad
export charge, xpos, ypos, pos
export VortexGroup, Dipole, Cluster
export Vortex, CoreShape, Ansatz, Exact
export detect_vortices_3d
export RecursiveClusterAlgorithm, RCAResult, recursive_cluster_algorithm

using .Core: rand_pointvortex

end
