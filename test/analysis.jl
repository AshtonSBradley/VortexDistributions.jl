abstract type Basis end

struct Oscillator <: Basis
    n::Int64
    ω::Float64
end

Oscillator(n::Int64) = Oscillator(n,1.0)
index(b::Oscillator) = 1:b.n
spectrum(b::Oscillator) = (0:b.n-1)*b.ω
qnumbers(b::Oscillator) = 0:b.n-1

function hermite(x,n)
    @assert n >= 1
    h1 = exp(-x^2 / 2) * (1 / π)^(1 / 4)
    n == 1 && (return h1)
    h2 = sqrt(2) * x * h1
    n == 2 && (return h2)
    h3 = 0.0
    for j ∈ 3:n
        h3 = sqrt(2/(j - 1)) * x * h2 - sqrt((j-2)/(j-1)) * h1
        h1 = h2; h2 = h3
    end
    n > 2 && (return h3)
end
hermite(x,n,ω) = hermite(√ω * x,n) * ω ^ ( 1 / 4 )

function (b::Oscillator)(x::Array{Float64,1})
    @unpack ω,n = b
    T = x * ones(1, n) |> zero
    T[:,1] = @. hermite(x,1,ω)
    n >= 2 && (T[:,2] = @. sqrt(2 * ω) * x * T[:,1])
    for j ∈ index(b)[3:end]
        m = qnumbers(b)[j]
        T[:,j] = @. sqrt(2 / m) * (√ω * x) * T[:,j-1] - sqrt((m-1) / m) * T[:,j-2]
    end
    return T
end

(b::Oscillator)(x::LinRange{Float64}) = b(x |> collect)
(b::Oscillator)(x::Float64,n::Int64) = hermite(x,n,b.ω)

function hermite_polar(x,n,ω)
    p = Progress(n,.1,"Precomputing...",20)   # minimum update interval: 1 second
    T = zeros(size(x)...,n)
    T[:,:,1] = @. hermite(x,1,ω)
    next!(p)
    n >= 2 && (T[:,:,2] = @. sqrt(2 * ω) * x * T[:,:,1])
    next!(p)
    for j in 3:n
        T[:,:,j] = @. sqrt(2 / (j-1)) * (√ω * x) * T[:,:,j-1] - sqrt((j-2) / (j-1)) * T[:,:,j-2]
    next!(p)
    end
    return T
end
hermite_polar(x,n) = hermite_polar(x,n,1.)

# filter cartesian data
function filterH(psi,x,N,ω)
    p = Progress(N,.1,"Filtering...",20)
    H = Oscillator(N,ω)(x)
    T = inv(H'*H)
    F̂ = T*H'*psi*H*T
    psiF = zero(psi)
        for m = 1:N
            next!(p)
            for n = 1:(N+1-m)
            @. psiF += H[:,m]*F̂[m,n]*H[:,n]'
            end
        end
    return psiF
end

# covert to polar coords
function slowpolar(psi,x,N,ω)
    H = Oscillator(N,ω)(x)
    T = inv(H'*H)
    F̂ = T*H'*psi*H*T

    Nx = length(x)
    θ = LinRange(0,2*pi,2*Nx)
    r = LinRange(0,last(x),Nx/2 |> Int)'

    psiP = zero(θ*r)
    for m = 1:N, n = 1:(N + 1 - m)
        @. psiP += hermite(r*cos(θ),m) * F̂[m,n] * hermite(r*sin(θ),n)
    end
    return psiP
end

function init_polar(x,N,ω)
    Nx = length(x)
    θ = LinRange(0,2*pi,2*Nx)
    r = LinRange(0,last(x),Nx/2 |> Int)'
    hx = hermite_polar((@. r*cos(θ)),N,ω)
    hy = hermite_polar((@. r*sin(θ)),N,ω)
    return hx,hy
end

function polar(psi,x,ω,Hx,Hy)
    N = size(Hx)[3]
    p = Progress(N,1,"Polar transform...",20)
    H = Oscillator(N,ω)(x)
    T = inv(H'*H)
    F̂ = T*H'*psi*H*T

    Nx = length(x)
    Nθ = 2*Nx
    Nr = Nx/2 |> Int
    θ = LinRange(0,2*pi,Nθ+1)[1:Nθ]
    r = LinRange(0,last(x),Nr)'

    psiFP = zeros(eltype(psi),Nθ,Nr)
    for m = 1:N
        next!(p)
        for n = 1:(N + 1 - m)
        hx = @view Hx[:,:,m]
        hy = @view Hy[:,:,n]
        @. psiFP += hx * F̂[m,n] * hy
    end
    end
    return psiFP
end
