using Test
using VortexDistributions

const RCA = VortexDistributions.RecursiveClusterAlgorithm

rca_pv(x, y, q) = PointVortex(x, y, q)

@testset "Dipoles" begin
    empty_result = RCA.get_dipoles(PointVortex[], 1.0, 1.0)
    @test isempty(empty_result.dipoles)
    @test isempty(empty_result.pairs)
    @test isempty(empty_result.mask)

    vortices = [rca_pv(0.10, 0.00, 1), rca_pv(0.12, 0.00, -1), rca_pv(0.80, 0.80, 1)]
    result = RCA.get_dipoles(vortices, 1.0, 1.0)
    @test result.pairs == [(1, 2)]
    @test length(result.dipoles) == 1
    @test result.mask == BitVector([true, true, false])

    periodic = [rca_pv(0.49, 0.00, 1), rca_pv(-0.49, 0.00, -1), rca_pv(0.00, 0.30, 1)]
    periodic_result = RCA.get_dipoles(periodic, 1.0, 1.0)
    @test periodic_result.pairs == [(1, 2)]
end

@testset "Seeds" begin
    vortices = [
        rca_pv(-0.30, 0.00, 1),
        rca_pv(-0.25, 0.00, 1),
        rca_pv(0.30, 0.00, -1),
        rca_pv(0.35, 0.00, -1),
    ]
    result = RCA.seed_clusters(vortices, 2.0, 2.0)
    @test result.pairs == [(1, 2), (3, 4)]
    @test result.mask == BitVector([true, true, true, true])
end

@testset "Cluster Growth" begin
    plus_vortices = [rca_pv(0.00, 0.00, 1), rca_pv(0.05, 0.00, 1), rca_pv(0.10, 0.00, 1), rca_pv(0.80, 0.00, -1)]
    plus_seeds = RCA.seed_clusters(plus_vortices, 2.0, 2.0)
    plus_clusters = RCA.grow_plus_clusters(plus_vortices, plus_seeds.pairs, 2.0, 2.0)
    @test length(plus_clusters) == 1
    @test length(only(plus_clusters).vortices) == 3

    minus_vortices = [rca_pv(-0.80, 0.00, 1), rca_pv(0.00, 0.00, -1), rca_pv(0.05, 0.00, -1), rca_pv(0.10, 0.00, -1)]
    minus_seeds = RCA.seed_clusters(minus_vortices, 2.0, 2.0)
    minus_clusters = RCA.grow_minus_clusters(minus_vortices, minus_seeds.pairs, 2.0, 2.0)
    @test length(minus_clusters) == 1
    @test length(only(minus_clusters).vortices) == 3
end

@testset "Driver" begin
    vortices = [
        rca_pv(-0.45, 0.00, 1),
        rca_pv(-0.43, 0.00, -1),
        rca_pv(0.00, 0.00, 1),
        rca_pv(0.05, 0.00, 1),
        rca_pv(0.10, 0.00, 1),
        rca_pv(0.70, 0.00, -1),
        rca_pv(0.75, 0.00, -1),
        rca_pv(0.40, 0.00, 1),
    ]

    result = RCA.recursive_cluster_algorithm(vortices, 2.0, 2.0)
    @test result isa RCA.RCAResult
    @test length(result.dipoles) == 1
    @test length(result.positive_clusters) == 1
    @test length(only(result.positive_clusters).vortices) == 3
    @test length(result.negative_clusters) == 1
    @test length(only(result.negative_clusters).vortices) == 2
    @test length(result.free_vortices) == 1
    @test only(result.free_vortices).xv == 0.40

    array_result = RCA.recursive_cluster_algorithm(
        [v.xv for v in vortices],
        [v.yv for v in vortices],
        [v.qv for v in vortices],
        2.0,
        2.0,
    )
    @test length(array_result.dipoles) == 1
    @test length(array_result.positive_clusters) == 1
end

@testset "Driver edge cases" begin
    dipole_only = [
        rca_pv(-0.10, 0.00, 1),
        rca_pv(-0.08, 0.00, -1),
    ]

    dipole_result = RCA.recursive_cluster_algorithm(dipole_only, 1.0, 1.0)
    @test length(dipole_result.dipoles) == 1
    @test isempty(dipole_result.positive_clusters)
    @test isempty(dipole_result.negative_clusters)
    @test isempty(dipole_result.free_vortices)

    array_dipole_result = RCA.recursive_cluster_algorithm(
        [v.xv for v in dipole_only],
        [v.yv for v in dipole_only],
        [v.qv for v in dipole_only],
        1.0,
        1.0,
    )
    @test length(array_dipole_result.dipoles) == 1
    @test isempty(array_dipole_result.free_vortices)

    @test_throws AssertionError RCA.recursive_cluster_algorithm([0.0], [0.0, 1.0], [1], 1.0, 1.0)
end
