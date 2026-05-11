"""
    DiskWindow(radius[, center])

Circular physical observation window for 2D point-vortex statistics.
"""
struct DiskWindow{T<:Real}
    radius::T
    center::Tuple{T,T}
end

function DiskWindow(radius::T) where {T<:Real}
    radius > 0 || throw(ArgumentError("disk radius must be positive"))
    return DiskWindow{T}(radius, (zero(T), zero(T)))
end

function DiskWindow(radius::Real, center::Tuple{<:Real,<:Real})
    return DiskWindow(promote(radius, center...)...)
end

function DiskWindow(radius::T, x0::T, y0::T) where {T<:Real}
    radius > 0 || throw(ArgumentError("disk radius must be positive"))
    return DiskWindow{T}(radius, (x0, y0))
end

"""
    RipleyKResult

Scale-resolved same-sign Ripley statistics. `K` is Ripley's K function,
`L = sqrt(K / π)` is Besag's L function, and `H = L - r` is positive
when same-sign vortices cluster relative to complete spatial randomness.
"""
struct RipleyKResult{T<:Real}
    r::Vector{T}
    K::Vector{T}
    L::Vector{T}
    H::Vector{T}
    charge::Int
    window::DiskWindow{T}
end

function _inside(window::DiskWindow, x::Real, y::Real)
    (; radius, center) = window
    x0, y0 = center
    return hypot(x - x0, y - y0) <= radius
end

_vortex_charge(vortex::PointVortex) = charge(vortex)

function _circle_fraction_inside_disk(window::DiskWindow, x::Real, y::Real, r::Real)
    (; radius, center) = window
    x0, y0 = center
    ρ = hypot(x - x0, y - y0)
    (!isfinite(ρ) || !isfinite(r)) && return zero(float(r))
    r <= radius - ρ && return one(float(r))
    (ρ == 0 || r == 0) && return zero(float(r))
    denominator = 2ρ * r
    (denominator == 0 || !isfinite(denominator)) && return zero(float(r))
    c = (radius^2 - ρ^2 - r^2) / denominator
    !isfinite(c) && return zero(float(c))
    c >= 1 && return one(float(c))
    c <= -1 && return zero(float(c))
    return 1 - acos(clamp(c, -1, 1)) / π
end

function _same_charge_positions(
    vortices::AbstractVector{PointVortex},
    q::Integer,
    window::DiskWindow,
)
    T = promote_type(Float64, typeof(window.radius))
    xs = T[]
    ys = T[]
    for vortex in vortices
        x = T(xpos(vortex))
        y = T(ypos(vortex))
        if charge(vortex) == q && _inside(window, x, y)
            push!(xs, x)
            push!(ys, y)
        end
    end
    return xs, ys
end

_rand(rng::Nothing) = rand()
_rand(rng) = rand(rng)
_randn(rng::Nothing) = randn()
_randn(rng) = randn(rng)

function _random_disk_vortices(rng, n::Integer, q::Integer, window::DiskWindow)
    (; radius, center) = window
    x0, y0 = center
    vortices = PointVortex[]
    sizehint!(vortices, n)
    for _ = 1:n
        ρ = radius * sqrt(_rand(rng))
        θ = 2π * _rand(rng)
        x = x0 + ρ * cos(θ)
        y = y0 + ρ * sin(θ)
        push!(vortices, PointVortex(x, y, q))
    end
    return vortices
end

"""
    ripley_k(vortices, r, window; charge = 1)

Compute the same-sign Ripley K function for `PointVortex` positions inside
`window`. Distances are physical Euclidean distances, and disk windows use an
isotropic boundary correction.
"""
function ripley_k(
    vortices::AbstractVector{PointVortex},
    r::AbstractVector{<:Real},
    window::DiskWindow;
    charge::Integer = 1,
)
    q = Int(charge)
    all(>=(0), r) || throw(ArgumentError("Ripley radii must be non-negative"))
    xs, ys = _same_charge_positions(vortices, q, window)
    T = promote_type(eltype(xs), eltype(r), typeof(window.radius))
    rs = T.(r)
    K = zeros(T, length(rs))
    n = length(xs)

    if n < 2
        L = zeros(T, length(rs))
        H = L .- rs
        return RipleyKResult(
            rs,
            K,
            L,
            H,
            q,
            DiskWindow{T}(T(window.radius), T.(window.center)),
        )
    end

    area = T(π) * T(window.radius)^2
    for i in eachindex(xs)
        xi = xs[i]
        yi = ys[i]
        for j in eachindex(xs)
            i == j && continue
            d = hypot(xi - xs[j], yi - ys[j])
            fraction = _circle_fraction_inside_disk(window, xi, yi, d)
            fraction = T(fraction)
            (!isfinite(fraction) || fraction <= eps(T)) && continue
            weight = inv(fraction)
            !isfinite(weight) && continue
            for k in eachindex(rs)
                d <= rs[k] && (K[k] += weight)
            end
        end
    end

    scale = area / (n * (n - 1))
    K .*= scale
    L = sqrt.(K ./ T(π))
    H = L .- rs
    return RipleyKResult(rs, K, L, H, q, DiskWindow{T}(T(window.radius), T.(window.center)))
end

function ripley_k(
    vortices::AbstractVector{PointVortex},
    r::Real,
    window::DiskWindow;
    charge::Integer = 1,
)
    return only(ripley_k(vortices, [r], window; charge).K)
end

"""
    besag_l(vortices, r, window; charge = 1)

Compute Besag's `L(r) = sqrt(K(r) / π)` for the same-sign vortex subset.
"""
function besag_l(
    vortices::AbstractVector{PointVortex},
    r::AbstractVector{<:Real},
    window::DiskWindow;
    charge::Integer = 1,
)
    return ripley_k(vortices, r, window; charge).L
end

function besag_l(
    vortices::AbstractVector{PointVortex},
    r::Real,
    window::DiskWindow;
    charge::Integer = 1,
)
    return only(besag_l(vortices, [r], window; charge))
end

"""
    ripley_csr_baseline(n, r, window; charge = 1, samples = 200, rng = nothing)

Estimate the homogeneous complete-spatial-randomness baseline in the same disk
window from `samples` random point sets with `n` same-sign vortices.
"""
function ripley_csr_baseline(
    n::Integer,
    r::AbstractVector{<:Real},
    window::DiskWindow;
    charge::Integer = 1,
    samples::Integer = 200,
    rng = nothing,
)
    samples > 0 || throw(ArgumentError("number of CSR samples must be positive"))
    rs = float.(r)
    mean_K = zeros(eltype(rs), length(rs))
    mean_L = similar(mean_K)
    mean_H = similar(mean_K)

    for _ = 1:samples
        vortices = _random_disk_vortices(rng, n, charge, window)
        result = ripley_k(vortices, r, window; charge)
        mean_K .+= result.K
        mean_L .+= result.L
        mean_H .+= result.H
    end

    mean_K ./= samples
    mean_L ./= samples
    mean_H ./= samples
    return (;
        r = rs,
        K = mean_K,
        L = mean_L,
        H = mean_H,
        charge = Int(charge),
        window,
        samples,
    )
end

function _positive_integral(r::AbstractVector, y::AbstractVector)
    integral = zero(eltype(y))
    for i = 2:length(r)
        Δr = r[i] - r[i-1]
        h = max((y[i] + y[i-1]) / 2, zero(integral))
        isfinite(h) || continue
        integral += h * Δr
    end
    return integral
end

"""
    ripley_excess(vortices, r, window; charge = 1, samples = 200, rng = nothing)

Compute `H_excess(r) = H_observed(r) - mean(H_CSR(r))` using homogeneous
random point sets in the same disk window. The returned `maximum`, `radius`,
and `integral` are evaluated on `H_excess`.
"""
function ripley_excess(
    vortices::AbstractVector{PointVortex},
    r::AbstractVector{<:Real},
    window::DiskWindow;
    charge::Integer = 1,
    samples::Integer = 200,
    rng = nothing,
)
    q = Int(charge)
    result = ripley_k(vortices, r, window; charge = q)
    n = count(
        vortex ->
            _vortex_charge(vortex) == q && _inside(window, xpos(vortex), ypos(vortex)),
        vortices,
    )
    baseline = ripley_csr_baseline(n, r, window; charge = q, samples, rng)
    H = result.H .- baseline.H
    isempty(H) &&
        return (; result, baseline, H, maximum = NaN, radius = NaN, integral = NaN)
    finite_indices = findall(isfinite, H)
    isempty(finite_indices) &&
        return (; result, baseline, H, maximum = NaN, radius = NaN, integral = NaN)
    index = finite_indices[argmax(H[finite_indices])]
    return (;
        result,
        baseline,
        H,
        maximum = H[index],
        radius = result.r[index],
        integral = _positive_integral(result.r, H),
    )
end

"""
    ripley_clustering(vortices, r, window; charge = 1)

Return the full Ripley result plus simple same-sign clustering summaries:
the maximum `H = L - r`, the radius where it occurs, and the positive
trapezoidal integral of `H`.
"""
function ripley_clustering(
    vortices::AbstractVector{PointVortex},
    r::AbstractVector{<:Real},
    window::DiskWindow;
    charge::Integer = 1,
)
    result = ripley_k(vortices, r, window; charge)
    isempty(result.H) && return (; result, maximum = NaN, radius = NaN, integral = NaN)
    index = argmax(result.H)
    integral = _positive_integral(result.r, result.H)
    return (; result, maximum = result.H[index], radius = result.r[index], integral)
end
