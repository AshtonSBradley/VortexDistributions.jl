using Test
using VortexDistributions

mutable struct TestRNG
    state::UInt64
end

function Base.rand(rng::TestRNG)
    rng.state = 6364136223846793005 * rng.state + 1442695040888963407
    return Float64(rng.state >> 11) * inv(Float64(1 << 53))
end

function Base.rand(rng::TestRNG, values::AbstractVector)
    return values[clamp(1 + floor(Int, rand(rng) * length(values)), 1, length(values))]
end

function Base.randn(rng::TestRNG)
    u1 = max(rand(rng), eps())
    u2 = rand(rng)
    return sqrt(-2log(u1)) * cos(2π * u2)
end

function test_random_disk_point(rng::TestRNG, radius::Real)
    ρ = radius * sqrt(rand(rng))
    θ = 2π * rand(rng)
    return ρ * cos(θ), ρ * sin(θ)
end

function test_random_disk_vortices(rng::TestRNG, n::Integer, q::Integer, radius::Real)
    return [PointVortex(test_random_disk_point(rng, radius)..., q) for _ in 1:n]
end

function test_clustered_disk_vortices(rng::TestRNG, n::Integer, q::Integer, radius::Real, σ::Real)
    centers = [(-0.35radius, 0.15radius), (0.30radius, -0.20radius)]
    vortices = PointVortex[]
    while length(vortices) < n
        x0, y0 = rand(rng, centers)
        x = x0 + σ * randn(rng)
        y = y0 + σ * randn(rng)
        hypot(x, y) <= radius && push!(vortices, PointVortex(x, y, q))
    end
    return vortices
end

@testset "Disk window construction" begin
    window = DiskWindow(10.0)
    @test window.radius == 10.0
    @test window.center == (0.0, 0.0)
    @test DiskWindow(10, (1, 2)).center == (1, 2)
    @test_throws ArgumentError DiskWindow(0.0)
end

@testset "Same-sign Ripley K" begin
    window = DiskWindow(10.0)
    vortices = [
        PointVortex(0.0, 0.0, 1),
        PointVortex(1.0, 0.0, 1),
        PointVortex(0.5, 0.0, -1),
        PointVortex(20.0, 0.0, 1),
    ]

    r = [0.5, 1.0, 2.0]
    result = ripley_k(vortices, r, window; charge = 1)

    @test result isa RipleyKResult
    @test result.charge == 1
    @test result.r == r
    @test result.K[1] == 0
    @test result.K[2] ≈ 100π
    @test result.K[3] ≈ 100π
    @test result.L[2] ≈ 10.0
    @test result.H[2] ≈ 9.0
    @test ripley_k(vortices, 1.0, window; charge = 1) ≈ 100π
    @test besag_l(vortices, 1.0, window; charge = 1) ≈ 10.0

    negative_result = ripley_k(vortices, r, window; charge = -1)
    @test all(iszero, negative_result.K)
    @test negative_result.H == -r
end

@testset "Disk edge correction" begin
    window = DiskWindow(1.0)
    vortices = [
        PointVortex(0.95, 0.0, 1),
        PointVortex(0.85, 0.0, 1),
    ]

    result = ripley_k(vortices, [0.1], window; charge = 1)
    @test only(result.K) > π
end

@testset "Ripley clustering summary" begin
    window = DiskWindow(10.0)
    vortices = [PointVortex(0.0, 0.0, 1), PointVortex(1.0, 0.0, 1)]
    summary = ripley_clustering(vortices, [0.5, 1.0, 2.0], window; charge = 1)

    @test summary.maximum ≈ 9.0
    @test summary.radius ≈ 1.0
    @test summary.integral > 0
    @test_throws ArgumentError ripley_k(vortices, [-1.0], window)
end

@testset "CSR-normalized Ripley excess" begin
    window = DiskWindow(1.0)
    r = collect(range(0.05, 0.5; length = 8))
    vortices = [PointVortex(0.0, 0.0, 1), PointVortex(0.05, 0.0, 1), PointVortex(0.1, 0.0, 1)]
    baseline = ripley_csr_baseline(3, r, window; samples = 8, rng = TestRNG(11))
    excess = ripley_excess(vortices, r, window; samples = 8, rng = TestRNG(11))

    @test length(baseline.H) == length(r)
    @test excess.H == excess.result.H .- excess.baseline.H
    @test excess.maximum == maximum(excess.H)
    @test excess.integral >= 0
end

@testset "Synthetic cluster peak regression" begin
    window = DiskWindow(1.0)
    σ = 0.07
    r = collect(range(0.02, 0.65; length = 64))
    data_rng = TestRNG(42)
    baseline_seed = 1

    homogeneous = test_random_disk_vortices(data_rng, 40, 1, window.radius)
    clustered = test_clustered_disk_vortices(data_rng, 40, 1, window.radius, σ)

    homogeneous_summary = ripley_excess(homogeneous, r, window; samples = 150, rng = TestRNG(baseline_seed))
    clustered_summary = ripley_excess(clustered, r, window; samples = 150, rng = TestRNG(baseline_seed))

    @test homogeneous_summary.maximum < 0.08
    @test clustered_summary.maximum > 0.35
    @test clustered_summary.integral > 20homogeneous_summary.integral
    @test clustered_summary.radius ≈ 2.86σ atol = 0.025
end
