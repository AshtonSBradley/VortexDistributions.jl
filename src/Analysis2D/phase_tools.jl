function zoom_grid(psi, x, y, xv, yv; win = 1, periodic = false)
    dx, dy = Δ(x), Δ(y)
    nx, ny = size(psi)
    ix = isapprox.(x, xv, atol = dx) |> findfirst
    iy = isapprox.(y, yv, atol = dy) |> findfirst
    ixw = (ix - win):(ix + win + 1)
    iyw = (iy - win):(iy + win + 1)
    if !periodic
        psiz, xz, yz = psi[ixw, iyw], x[ixw], y[iyw]
    else
        ixp, iyp = mod1.(ixw, nx), mod1.(iyw, ny)
        psiz = psi[ixp, iyp]
        xz = (x[ix] - win * dx):dx:(x[ix] + (win + 1) * dx)
        yz = (y[iy] - win * dy):dy:(y[iy] + (win + 1) * dy)
    end
    return psiz, xz, yz
end

function zoom_interp(psi, x, y, xv, yv; win = 1, nz = 30, periodic = false)
    psiz, xz, yz = zoom_grid(psi, x, y, xv, yv, win = win, periodic = periodic)
    psi_itp = interpolate((xz, yz), psiz, Gridded(Linear()))
    xint, yint = LinRange(xz[1], xz[end], nz), LinRange(yz[1], yz[end], nz)
    psi_int = psi_itp(xint, yint)
    return psi_int, xint, yint
end

circ_mask(x, y, R) = hypot(x, y) <= R

function keep_vortices(vort, mask = (x, y) -> circ_mask(x, y, 1.0))
    vortm = PointVortex[]
    for v in vort
        x, y = pos(v)
        if mask(x, y)
            push!(vortm, v)
        end
    end
    return vortm
end

function phase_jumps(phase, dim = 1)
    @assert (dim == 1 || dim == 2)
    s1, s2 = size(phase)
    pdiff = zero(phase)

    if dim == 1
        for j in 1:s2
            @inbounds dϕ = phase[1, j] - phase[s1, j]
            @inbounds abs(dϕ) > π && (pdiff[1, j] += sign(dϕ))
            for i in 2:s1
                @inbounds dϕ = phase[i, j] - phase[i - 1, j]
                @inbounds abs(dϕ) > π && (pdiff[i, j] += sign(dϕ))
            end
        end
    else
        for j in 2:s2
            for i in 1:s1
                j == 2 && (@inbounds dϕ = phase[i, 1] - phase[i, s2];
                    @inbounds abs(dϕ) > π && (pdiff[i, 1] += sign(dϕ)))
                @inbounds dϕ = phase[i, j] - phase[i, j - 1]
                @inbounds abs(dϕ) > π && (pdiff[i, j] += sign(dϕ))
            end
        end
    end
    return pdiff
end

function phase_jumps!(pdiff, phase, dim = 1)
    @assert (dim == 1 || dim == 2)
    s1, s2 = size(phase)

    if dim == 1
        for j in 1:s2
            @inbounds dϕ = phase[1, j] - phase[s1, j]
            @inbounds abs(dϕ) > π && (pdiff[1, j] += sign(dϕ))
            for i in 2:s1
                @inbounds dϕ = phase[i, j] - phase[i - 1, j]
                @inbounds abs(dϕ) > π && (pdiff[i, j] += sign(dϕ))
            end
        end
    else
        for j in 2:s2
            for i in 1:s1
                j == 2 && (@inbounds dϕ = phase[i, 1] - phase[i, s2];
                    @inbounds abs(dϕ) > π && (pdiff[i, 1] += sign(dϕ)))
                @inbounds dϕ = phase[i, j] - phase[i, j - 1]
                @inbounds abs(dϕ) > π && (pdiff[i, j] += sign(dϕ))
            end
        end
    end
end

function unwrap(phase::Array{Float64, 2}, dim = 1)
    @assert (dim == 1 || dim == 2)
    s1, s2 = size(phase)
    uphase = copy(phase)

    if dim == 1
        for j in 1:s2
            @inbounds (k = uphase[1, j] - uphase[s1, j]; k > π) && (uphase[1, j] -= 2π * div(k, 2pi, RoundNearest))
            @inbounds (k = uphase[1, j] - uphase[s1, j]; k < -π) && (uphase[1, j] += 2π * div(-k, 2pi, RoundNearest))
            for i in 2:s1
                @inbounds (k = uphase[i, j] - uphase[i - 1, j]; k > π) && (uphase[i, j] -= 2π * div(k, 2pi, RoundNearest))
                @inbounds (k = uphase[i, j] - uphase[i - 1, j]; k < -π) && (uphase[i, j] += 2π * div(-k, 2pi, RoundNearest))
            end
        end
    else
        for j in 2:s2
            for i in 1:s1
                j == 2 && (@inbounds (k = uphase[i, 1] - uphase[i, s2]; k > π) && (uphase[i, 1] -= 2π * div(k, 2pi, RoundNearest));
                    @inbounds (k = uphase[i, 1] - uphase[i, s2]; k < -π) && (uphase[i, 1] += 2π * div(-k, 2pi, RoundNearest)))
                @inbounds (k = uphase[i, j] - uphase[i, j - 1]; k > π) && (uphase[i, j] -= 2π * div(k, 2pi, RoundNearest))
                @inbounds (k = uphase[i, j] - uphase[i, j - 1]; k < -π) && (uphase[i, j] += 2π * div(-k, 2pi, RoundNearest))
            end
        end
    end

    return uphase
end

function unwrap!(uphase::Array{Float64, 2}, phase::Array{Float64, 2}, dim = 1)
    @assert (dim == 1 || dim == 2)
    s1, s2 = size(phase)
    uphase .= phase

    if dim == 1
        for j in 1:s2
            @inbounds (k = uphase[1, j] - uphase[s1, j]; k > π) && (uphase[1, j] -= 2π * div(k, 2pi, RoundNearest))
            @inbounds (k = uphase[1, j] - uphase[s1, j]; k < -π) && (uphase[1, j] += 2π * div(-k, 2pi, RoundNearest))
            for i in 2:s1
                @inbounds (k = uphase[i, j] - uphase[i - 1, j]; k > π) && (uphase[i, j] -= 2π * div(k, 2pi, RoundNearest))
                @inbounds (k = uphase[i, j] - uphase[i - 1, j]; k < -π) && (uphase[i, j] += 2π * div(-k, 2pi, RoundNearest))
            end
        end
    else
        for j in 2:s2
            for i in 1:s1
                j == 2 && (@inbounds (k = uphase[i, 1] - uphase[i, s2]; k > π) && (uphase[i, 1] -= 2π * div(k, 2pi, RoundNearest));
                    @inbounds (k = uphase[i, 1] - uphase[i, s2]; k < -π) && (uphase[i, 1] += 2π * div(-k, 2pi, RoundNearest)))
                @inbounds (k = uphase[i, j] - uphase[i, j - 1]; k > π) && (uphase[i, j] -= 2π * div(k, 2pi, RoundNearest))
                @inbounds (k = uphase[i, j] - uphase[i, j - 1]; k < -π) && (uphase[i, j] += 2π * div(-k, 2pi, RoundNearest))
            end
        end
    end
end

function unwrap(phase::Array{Float64, 1})
    uphase = copy(phase)
    s1 = length(phase)
    @inbounds (k = uphase[1] - uphase[s1]; k > π) && (uphase[1] -= 2π * div(k, 2pi, RoundNearest))
    @inbounds (k = uphase[1] - uphase[s1]; k < -π) && (uphase[1] += 2π * div(-k, 2pi, RoundNearest))
    for i in 2:s1
        @inbounds (k = uphase[i] - uphase[i - 1]; k > π) && (uphase[i] -= 2π * div(k, 2pi, RoundNearest))
        @inbounds (k = uphase[i] - uphase[i - 1]; k < -π) && (uphase[i] += 2π * div(-k, 2pi, RoundNearest))
    end
    return uphase
end

function unwrap!(uphase::Array{Float64, 1}, phase::Array{Float64, 1})
    uphase .= phase
    s1 = length(phase)
    @inbounds (k = uphase[1] - uphase[s1]; k > π) && (uphase[1] -= 2π * div(k, 2pi, RoundNearest))
    @inbounds (k = uphase[1] - uphase[s1]; k < -π) && (uphase[1] += 2π * div(-k, 2pi, RoundNearest))
    for i in 2:s1
        @inbounds (k = uphase[i] - uphase[i - 1]; k > π) && (uphase[i] -= 2π * div(k, 2pi, RoundNearest))
        @inbounds (k = uphase[i] - uphase[i - 1]; k < -π) && (uphase[i] += 2π * div(-k, 2pi, RoundNearest))
    end
end
