using Printf
using Random
using VortexDistributions

function clustered_disk_vortices(
    rng::AbstractRNG,
    n::Integer,
    q::Integer,
    radius::Real,
    cluster_width::Real;
    centers = [(-0.35radius, 0.15radius), (0.30radius, -0.20radius)],
)
    vortices = PointVortex[]
    while length(vortices) < n
        x0, y0 = rand(rng, centers)
        x = x0 + cluster_width * randn(rng)
        y = y0 + cluster_width * randn(rng)
        hypot(x, y) <= radius && push!(vortices, PointVortex(x, y, q))
    end
    return vortices
end

function summarize(values)
    sorted = sort(values)
    n = length(sorted)
    median = isodd(n) ? sorted[(n+1)÷2] : (sorted[n÷2] + sorted[n÷2+1]) / 2
    lo = sorted[max(1, round(Int, 0.16n))]
    hi = sorted[min(n, round(Int, 0.84n))]
    return median, lo, hi
end

function calibrate_peak(
    rng::AbstractRNG;
    radius = 1.0,
    n = 40,
    widths = [0.04, 0.06, 0.08, 0.10, 0.12],
    realizations = 80,
    baseline_samples = 200,
)
    window = DiskWindow(radius)
    r = collect(range(0.02radius, 0.65radius; length = 96))

    println(
        "Ripley-excess peak calibration for two Gaussian same-sign clusters in a unit disk",
    )
    println("Each cluster has coordinate width σ; table reports median and 16-84% range.")
    println()
    println(
        "σ       median r_peak   16-84% r_peak       median r_peak/σ   median r_peak/(2σ)",
    )

    for σ in widths
        peaks = Float64[]
        for _ = 1:realizations
            vortices = clustered_disk_vortices(rng, n, 1, radius, σ)
            summary = ripley_excess(
                vortices,
                r,
                window;
                charge = 1,
                samples = baseline_samples,
                rng,
            )
            push!(peaks, summary.radius)
        end

        median, lo, hi = summarize(peaks)
        @printf(
            "%0.3f   %0.4f          [%0.4f, %0.4f]      %0.3f             %0.3f\n",
            σ,
            median,
            lo,
            hi,
            median / σ,
            median / (2σ),
        )
    end
end

calibrate_peak(MersenneTwister(7))
