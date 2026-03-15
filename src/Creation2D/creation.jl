scalar_ansatz(x) = sqrt(x^2 / (1 + x^2))

Ansatz() = Ansatz(ψa, 1.0, Λ)

function (core::Ansatz)(x)
    @unpack f, ξ, Λ = core
    return f(Λ * x / ξ)
end

(core::Ansatz)(x, y) = core(hypot(x, y))

Exact(ξ::Float64) = Exact(ψi, ξ)
Exact() = Exact(1.0)

function (core::Exact)(x)
    @unpack ξ, f = core
    return f(x / ξ)
end

(core::Exact)(x, y) = core(hypot(x, y))

ScalarVortex(vort::PointVortex) = ScalarVortex(Exact(), vort)

function (s::ScalarVortex{T})(x, y) where {T <: CoreShape}
    @unpack xv, yv, qv = s.vort
    return s.core(x - xv, y - yv) * exp(im * qv * atan(y - yv, x - xv))
end

ScalarVortex(ξ::Float64, pv::Array{PointVortex, 1}) = ScalarVortex.([Exact(ξ)], pv::Array{PointVortex, 1})
ScalarVortex(ξ::Array{Float64, 1}, pv::Array{PointVortex, 1}) = @. ScalarVortex(Exact(ξ), pv::Array{PointVortex, 1})
ScalarVortex(ξ::Float64, pv::PointVortex) = ScalarVortex(ξ, [pv])
ScalarVortex(pv::Array{PointVortex, 1}) = ScalarVortex(1.0, pv)

ChannelVortex(vort::PointVortex) = ChannelVortex(Exact(), vort)

function (s::ChannelVortex{T})(x, y, z::Real) where {T <: CoreShape}
    @unpack xv, yv, qv = s.vort
    return s.core(x - xv, y - yv) * exp(im * qv * (real(1im * log((exp(π * (x + 1im * (y + z / 2)) / z)
        - exp(π * (xv - 1im * (yv + z / 2)) / z)) / (exp(π * (x + 1im * (y + z / 2)) / z)
        - exp(π * (xv + 1im * (yv + z / 2)) / z) + 1e-10)))))
end

rand_scalarvortex() = ScalarVortex(rand_pointvortex())
rand_scalarvortex(n) = ScalarVortex.(rand_pointvortex(n))
rand_scalarvortex(psi::Field) = ScalarVortex(Exact(), rand_pointvortex(psi))
rand_scalarvortex(n, psi::Field) = ScalarVortex.([Exact()], rand_pointvortex(n, psi))

rand_vortex() = rand_scalarvortex()
rand_vortex(n) = rand_scalarvortex(n)
rand_vortex(n, psi::Field) = rand_scalarvortex(n, psi)

function vortex!(psi::F, vort::ScalarVortex{T}) where {T <: CoreShape, F <: Field}
    @unpack ψ, x, y = psi
    ψ .*= vort.(x, y')
    @pack! psi = ψ
end

function vortex!(psi::F, vort::S) where {F <: Field, S <: Vortex}
    vortex!(psi, ScalarVortex(vort))
end

function vortex!(psi::F, vort::Array{S}) where {F <: Field, S <: Vortex}
    for j in eachindex(vort)
        vortex!(psi, vort[j])
    end
end

function vortex!(psi::F, vort::ChannelVortex{T}, D) where {T <: CoreShape, F <: Field}
    @unpack ψ, x, y = psi
    @assert typeof(D) <: Real "Channel length cannot be imaginary"
    ψ .*= vort.(x, y', D)
    @pack! psi = ψ
end

function rand_vortexfield(n)
    Nx = 400
    Ny = 400
    Lx = 200
    Ly = 200
    x = LinRange(-Lx / 2, Ly / 2, Nx)
    y = LinRange(-Ly / 2, Ly / 2, Ny)

    psi0 = one.(x * y') |> complex
    psi = Torus(psi0, x, y)
    vort = rand_vortex(n, psi)
    vortex!(psi, vort)
    return psi, PointVortex(vort)
end

H(x) = x >= 0.0 ? 1.0 : 0.0
shift(x, y) = x - y
tans(x, xk) = tan((shift(x, xk) - π) * 0.5)
tanhs(x, xk, j) = tanh((shift(x, xk) + 2 * π * j) * 0.5)

function kernel(x, y, xp, yp, xn, yn, j)
    return atan(tanhs(y, yn, j) * tans(x, xn)) -
        atan(tanhs(y, yp, j) * tans(x, xp))
end

function thetad(x, y, xp, yp, xn, yn)
    s = 0.0
    for j in -5:5
        s += kernel(x, y, xp, yp, xn, yn, j)
    end
    s += π * (H(shift(x, xp)) - H(shift(x, xn))) - y * (xp - xn) / (2 * π)
    return s - x * H(abs(yp - yn) - π) / (2 * pi) + y * H(abs(xp - xn) - π) / (2 * pi)
end

function dipole_phase(x, y, xp, yp, xn, yn)
    Lx = x[end] - x[1]
    Ly = y[end] - y[1]
    return @. thetad(x * 2 * pi / Lx, y' * 2 * pi / Ly, xp * 2 * pi / Lx, yp * 2 * pi / Ly, xn * 2 * pi / Lx, yn * 2 * pi / Ly)
end

function periodic_dipole!(psi::F, dip::Vector{ScalarVortex{T}}) where {T <: CoreShape, F <: Field}
    @assert length(dip) == 2
    @assert dip[1].vort.qv + dip[2].vort.qv == 0
    @unpack ψ, x, y = psi
    (dip[1].vort.qv > 0) ? (jp = 1; jn = 2) : (jp = 2; jn = 1)
    vp = vortex_array(dip[jp].vort)[1:2]
    vn = vortex_array(dip[jn].vort)[1:2]
    @. ψ *= abs(dip[jn](x, y') * dip[jp](x, y'))
    ψ .*= exp.(im * dipole_phase(x, y, vp..., vn...))
    @pack! psi = ψ
end

function gpecore_exact(K, L = 2, N = 100, R = K)
    z, Dz = getChebDMatrix(N)
    _, D2z = getChebD2Matrix(N)
    z = vec(reverse(z, dims = 2))
    z = z * L ./ 2
    Dz = -2 * Dz ./ L
    D2z = 4 * D2z ./ L.^2

    y = @. R * (1 + z) / (1 - z)

    ψ = @. y / hypot(1 / Λ, y)
    ψ[1] = 0
    ψ[end] = 1

    Q = vec(z.^2 .- 2 * z .+ 1)
    Qmat = repeat(Q, 1, N)
    Zmat = repeat(2 * (z .- 1), 1, N)
    Ymat = repeat(1 ./ y, 1, N)

    residuals = -0.5 * ((Q .* (D2z * ψ) + 2 * (z .- 1) .* (Dz * ψ)) .* Q / (4 * R^2)
        + (Q / (2 * R) ./ y) .* (Dz * ψ)) + 0.5 * K^2 * ψ ./ y.^2 + ψ.^3 .- ψ
    residuals[1] = 0
    residuals[end] = 0

    while sum(abs.(residuals).^2) > 1e-12
        residuals = -0.5 * ((Q .* (D2z * ψ) + 2 * (z .- 1) .* (Dz * ψ)) .* Q / (4 * R^2)
            + (Q / (2 * R) ./ y) .* (Dz * ψ)) + 0.5 * K^2 * ψ ./ y.^2 + ψ.^3 .- ψ
        residuals[1] = 0
        residuals[end] = 0

        Jacobi = -0.5 * ((Qmat .* D2z + Zmat .* Dz) .* (Qmat / (4 * R^2))
            + (Qmat / (2 * R) .* Ymat) .* Dz) + diagm(0 => 0.5 * K^2 ./ y.^2) + diagm(0 => 3 * ψ.^2 .- 1)

        Jacobi[1, :] = [1 zeros(1, N - 1)]
        Jacobi[N, :] = [zeros(1, N - 1) 1]

        Δ = vec(-4 / 7 * (Jacobi \ residuals))

        ψ = ψ + Δ
        ψ[1] = 0
        ψ[end] = 1
    end
    res = norm(residuals)^2
    return y, ψ, res
end

function chebdif(N, M)
    @assert 0 < M <= N - 1

    n1, n2 = (floor(N / 2) |> Int, ceil(N / 2) |> Int)

    k = collect(0.0:N - 1)'
    th = k * π / (N - 1)

    x = sin.(π * collect(N - 1:-2:1 - N)' ./ (2 * (N - 1)))

    T = repeat(th / 2, N, 1)'
    DX = 2 * sin.(T' .+ T) .* sin.(T' .- T)
    DX = [DX[1:n1, :]; -reverse(reverse(DX[1:n2, :], dims = 2), dims = 1)]
    DX[diagind(DX)] .= 1

    C = [(-1.0)^(i - j) for i in 0:(N - 1), j in 0:(N - 1)]
    C[1, :] *= 2
    C[N, :] *= 2
    C[:, 1] /= 2
    C[:, N] /= 2

    Z = 1 ./ DX
    Z[diagind(Z)] .= 0

    D = Matrix{Float64}(I, N, N)
    DM = zeros(N, N, M)

    for ell in 1:M
        D = ell * Z .* (C .* repeat(diag(D), 1, N) - D)
        D[diagind(D)] .= -sum(D, dims = 2)
        DM[:, :, ell] = D
    end
    return x, DM
end

function getChebDMatrix(N)
    z, DM = chebdif(N, 1)
    return z, DM[:, :, 1]
end

function getChebD2Matrix(N)
    z, DM = chebdif(N, 2)
    return z, DM[:, :, 2]
end
