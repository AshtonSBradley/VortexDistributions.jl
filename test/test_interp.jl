using Test
using VortexDistributions

N = 32
L = 8.0
dx = L / N
x = collect(range(-L / 2, step = dx, length = N))
y = copy(x)
psi = Torus(ones(ComplexF64, N, N), x, y)

grid_patch, grid_x, grid_y = zoom_grid(psi.ψ, psi.x, psi.y, [0.0], [0.0]; win = 1, periodic = false)
interp_patch, interp_x, interp_y = zoom_interp(psi.ψ, psi.x, psi.y, [0.0], [0.0]; win = 1, nz = 10, periodic = false)

@test size(grid_patch) == (4, 4)
@test all(grid_patch .== 1.0 + 0.0im)
@test grid_x == [-0.5, -0.25, 0.0, 0.25]
@test grid_y == [-0.5, -0.25, 0.0, 0.25]

@test size(interp_patch) == (10, 10)
@test all(interp_patch .≈ (1.0 + 0.0im))
@test length(interp_x) == 10
@test length(interp_y) == 10
