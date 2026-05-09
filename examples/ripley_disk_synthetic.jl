using Printf
using Random
using VortexDistributions

function random_disk_point(rng::AbstractRNG, radius::Real)
    ρ = radius * sqrt(rand(rng))
    θ = 2π * rand(rng)
    return ρ * cos(θ), ρ * sin(θ)
end

function random_disk_vortices(rng::AbstractRNG, n::Integer, q::Integer, radius::Real)
    return [PointVortex(random_disk_point(rng, radius)..., q) for _ = 1:n]
end

function clustered_disk_vortices(
    rng::AbstractRNG,
    n::Integer,
    q::Integer,
    radius::Real;
    centers = [(-0.35radius, 0.15radius), (0.30radius, -0.20radius)],
    cluster_width = 0.07radius,
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

function print_summary(name, summary)
    @printf(
        "%-12s max(H_excess) = %8.4f at r = %6.3f, positive integral = %8.4f\n",
        name,
        summary.maximum,
        summary.radius,
        summary.integral,
    )
end

function points_path(x, y, xmin, xmax, ymin, ymax, width, height, margin)
    xs = @. margin + (x - xmin) / (xmax - xmin) * (width - 2margin)
    ys = @. height - margin - (y - ymin) / (ymax - ymin) * (height - 2margin)
    return join(
        ("$(round(xs[i], digits = 2)),$(round(ys[i], digits = 2))" for i in eachindex(xs)),
        " ",
    )
end

function xcoord(x, xmin, xmax, width, margin)
    return margin + (x - xmin) / (xmax - xmin) * (width - 2margin)
end

function ycoord(y, ymin, ymax, height, margin)
    return height - margin - (y - ymin) / (ymax - ymin) * (height - 2margin)
end

function write_svg(path, r, homogeneous, clustered; calibrated_peak_scale)
    width = 780
    height = 460
    margin = 56
    ymin = min(minimum(homogeneous), minimum(clustered), -0.05)
    ymax = max(maximum(homogeneous), maximum(clustered), 0.05)
    pad = 0.08 * (ymax - ymin)
    ymin -= pad
    ymax += pad
    xmin, xmax = extrema(r)
    zero_y = ycoord(0, ymin, ymax, height, margin)
    hom_path = points_path(r, homogeneous, xmin, xmax, ymin, ymax, width, height, margin)
    clu_path = points_path(r, clustered, xmin, xmax, ymin, ymax, width, height, margin)
    scale_x = xcoord(calibrated_peak_scale, xmin, xmax, width, margin)
    xticks = collect(0.1:0.1:0.6)
    yticks = collect(0.0:0.1:0.5)

    open(path, "w") do io
        println(
            io,
            """<svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" viewBox="0 0 $width $height">""",
        )
        println(io, """<rect width="100%" height="100%" fill="white"/>""")
        println(
            io,
            """<text x="24" y="34" font-family="Helvetica,Arial,sans-serif" font-size="20" fill="#222">Same-sign Ripley excess in a disk</text>""",
        )
        println(
            io,
            """<line x1="$margin" y1="$(height - margin)" x2="$(width - margin)" y2="$(height - margin)" stroke="#333"/>""",
        )
        println(
            io,
            """<line x1="$margin" y1="$margin" x2="$margin" y2="$(height - margin)" stroke="#333"/>""",
        )
        println(
            io,
            """<line x1="$margin" y1="$zero_y" x2="$(width - margin)" y2="$zero_y" stroke="#999" stroke-dasharray="5 5"/>""",
        )
        for tick in xticks
            x = xcoord(tick, xmin, xmax, width, margin)
            println(
                io,
                """<line x1="$x" y1="$(height - margin)" x2="$x" y2="$(height - margin + 6)" stroke="#333"/>""",
            )
            println(
                io,
                """<text x="$x" y="$(height - margin + 24)" text-anchor="middle" font-family="Helvetica,Arial,sans-serif" font-size="12" fill="#222">$(round(tick, digits = 1))</text>""",
            )
        end
        for tick in yticks
            y = ycoord(tick, ymin, ymax, height, margin)
            println(
                io,
                """<line x1="$(margin - 6)" y1="$y" x2="$margin" y2="$y" stroke="#333"/>""",
            )
            println(
                io,
                """<text x="$(margin - 10)" y="$(y + 4)" text-anchor="end" font-family="Helvetica,Arial,sans-serif" font-size="12" fill="#222">$(round(tick, digits = 1))</text>""",
            )
        end
        println(
            io,
            """<line x1="$scale_x" y1="$margin" x2="$scale_x" y2="$(height - margin)" stroke="#111" stroke-width="2.5" stroke-dasharray="6 4"/>""",
        )
        println(io, """<circle cx="$scale_x" cy="$(height - margin)" r="4" fill="#111"/>""")
        println(
            io,
            """<polyline fill="none" stroke="#4b6cb7" stroke-width="3" points="$hom_path"/>""",
        )
        println(
            io,
            """<polyline fill="none" stroke="#c43c39" stroke-width="3" points="$clu_path"/>""",
        )
        println(
            io,
            """<text x="$(width - 240)" y="72" font-family="Helvetica,Arial,sans-serif" font-size="15" fill="#4b6cb7">homogeneous - CSR</text>""",
        )
        println(
            io,
            """<text x="$(width - 240)" y="96" font-family="Helvetica,Arial,sans-serif" font-size="15" fill="#c43c39">clustered - CSR</text>""",
        )
        println(
            io,
            """<text x="$scale_x" y="$(margin - 14)" text-anchor="middle" font-family="Helvetica,Arial,sans-serif" font-size="14" font-weight="bold" fill="#111">reference r = $(round(calibrated_peak_scale, digits = 3))</text>""",
        )
        println(
            io,
            """<text x="$(width / 2 - 40)" y="$(height - 16)" font-family="Helvetica,Arial,sans-serif" font-size="15" fill="#222">r</text>""",
        )
        println(
            io,
            """<text x="14" y="$(height / 2)" transform="rotate(-90 14,$(height / 2))" font-family="Helvetica,Arial,sans-serif" font-size="15" fill="#222">H_excess(r)</text>""",
        )
        println(io, "</svg>")
    end
end

rng = MersenneTwister(42)
window = DiskWindow(1.0)
r = collect(range(0.02, 0.65; length = 64))
cluster_width = 0.07window.radius

homogeneous = random_disk_vortices(rng, 40, 1, window.radius)
clustered = clustered_disk_vortices(rng, 40, 1, window.radius; cluster_width)

homogeneous_summary = ripley_excess(
    homogeneous,
    r,
    window;
    charge = 1,
    samples = 300,
    rng = MersenneTwister(1),
)
clustered_summary =
    ripley_excess(clustered, r, window; charge = 1, samples = 300, rng = MersenneTwister(1))

println("Same-sign Ripley clustering in a disk window")
print_summary("homogeneous", homogeneous_summary)
print_summary("clustered", clustered_summary)

plot_path = joinpath(@__DIR__, "ripley_disk_synthetic.svg")
write_svg(
    plot_path,
    r,
    homogeneous_summary.H,
    clustered_summary.H;
    calibrated_peak_scale = 2.86cluster_width,
)
println("Wrote figure: ", plot_path)

println("\nSample H_excess(r) = H_observed(r) - mean(H_CSR(r)) values")
println("r      homogeneous   clustered")
for index = 1:8:length(r)
    @printf(
        "%5.3f  %11.4f  %10.4f\n",
        r[index],
        homogeneous_summary.H[index],
        clustered_summary.H[index],
    )
end
