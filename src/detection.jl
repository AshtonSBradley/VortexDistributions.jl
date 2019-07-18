"""
    jumps = phasejumps(phase,dim)

Count phase jumps greater than π in `phase` along dimension `dim`
"""
function phasejumps(phase,dim=1)
    @assert (dim==1 || dim==2)
    Nx,Ny = size(phase)
    pdiff = zero(phase)

    if dim == 1
    @inbounds for j in 1:Ny
        @inbounds for i in 2:Nx
            Δ = phase[i,j] - phase[i-1,j]
            abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end
            Δ = phase[1,j] - phase[Nx,j]
            abs(Δ) > π && (pdiff[1,j] += sign(Δ))
    end

    elseif dim == 2
    @inbounds for j in 2:Ny
        @inbounds for i in 1:Nx
            Δ = phase[i,j] - phase[i,j-1]
            abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end
    end
        @inbounds for i in 1:Nx
            Δ = phase[i,1] - phase[i,Ny]
            abs(Δ) > π && (pdiff[i,1] += sign(Δ))
        end
    end
  return pdiff
end

function phasejumps!(pdiff,phase,dim=1)
    @assert (dim==1 || dim==2)
    Nx,Ny = size(phase)

    if dim == 1
    @inbounds for j in 1:Ny
        @inbounds for i in 2:Nx
            Δ = phase[i,j] - phase[i-1,j]
            abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end
            Δ = phase[1,j] - phase[Nx,j]
            abs(Δ) > π && (pdiff[1,j] += sign(Δ))
    end

    elseif dim == 2
    @inbounds for j in 2:Ny
        @inbounds for i in 1:Nx
            Δ = phase[i,j] - phase[i,j-1]
            abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end

    end
        @inbounds for i in 1:Nx
            Δ = phase[i,1] - phase[i,Ny]
            abs(Δ) > π && (pdiff[i,1] += sign(Δ))
        end
    end
end

"""
`unwrapped = unwrap(phase,dim=1)`

Unwraps 2d array `phase` along dimension `dim`, acting periodically to give back array of same size as `phase`.

`unwrap!(unwrapped,phase,dim)` writes in-place to `unwrapped`.
"""
function unwrap(phase::Array{Float64,2},dim=1)
    @assert (dim==1 || dim==2)
    Nx,Ny = size(phase)
    uphase = copy(phase)

    if dim == 1
    @inbounds for j in 1:Ny
        @inbounds for i in 2:Nx
        (uphase[i,j] - uphase[i-1,j] >= π) && (uphase[i,j] -= 2π)
        (uphase[i,j] - uphase[i-1,j] <= -π) && (uphase[i,j] += 2π)
        end
        (uphase[1,j] - uphase[Nx,j] >= π) && (uphase[1,j] -= 2π)
        (uphase[1,j] - uphase[Nx,j] <= -π) && (uphase[1,j] += 2π)
    end

    elseif dim == 2
    @inbounds for j in 2:Ny
        @inbounds for i in 1:Nx
        (uphase[i,j] - uphase[i,j-1] >= π) && (uphase[i,j] -= 2π)
        (uphase[i,j] - uphase[i,j-1] <= -π) && (uphase[i,j] += 2π)
        end

    end
        @inbounds for i in 1:Nx
        (uphase[i,1] - uphase[i,Ny] >= π) && (uphase[i,1] -= 2π)
        (uphase[i,1] - uphase[i,Ny] <= -π) && (uphase[i,1] += 2π)
        end
    end

    return uphase
end

function unwrap!(uphase::Array{Float64,2},phase::Array{Float64,2},dim=1)
    @assert (dim==1 || dim==2)
    Nx,Ny = size(phase)
    uphase .= phase

    if dim == 1
    @inbounds for j in 1:Ny
        @inbounds for i in 2:Nx
        (uphase[i,j] - uphase[i-1,j] >= π) && (uphase[i,j] -= 2π)
        (uphase[i,j] - uphase[i-1,j] <= -π) && (uphase[i,j] += 2π)
        end
        (uphase[1,j] - uphase[Nx,j] >= π) && (uphase[1,j] -= 2π)
        (uphase[1,j] - uphase[Nx,j] <= -π) && (uphase[1,j] += 2π)
    end

    elseif dim == 2
    @inbounds for j in 2:Ny
        @inbounds for i in 1:Nx
        (uphase[i,j] - uphase[i,j-1] >= π) && (uphase[i,j] -= 2π)
        (uphase[i,j] - uphase[i,j-1] <= -π) && (uphase[i,j] += 2π)
        end

    end
        @inbounds for i in 1:Nx
        (uphase[i,1] - uphase[i,Ny] >= π) && (uphase[i,1] -= 2π)
        (uphase[i,1] - uphase[i,Ny] <= -π) && (uphase[i,1] += 2π)
        end
    end
end

function unwrap(phase::Array{Float64,1})
    uphase = copy(phase)
    Nx = length(phase)
        @inbounds for i in 2:Nx
        (uphase[i] - uphase[i-1] >= π) && (uphase[i] -= 2π)
        (uphase[i] - uphase[i-1] <= -π) && (uphase[i] += 2π)
        end
        (uphase[1] - uphase[Nx] >= π) && (uphase[1] -= 2π)
        (uphase[1] - uphase[Nx] <= -π) && (uphase[1] += 2π)
        return uphase
end

function unwrap!(uphase::Array{Float64,1},phase::Array{Float64,1})
    Nx = length(phase)
        @inbounds for i in 2:Nx
        (uphase[i] - uphase[i-1] >= π) && (uphase[i] -= 2π)
        (uphase[i] - uphase[i-1] <= -π) && (uphase[i] += 2π)
        end
        (uphase[1] - uphase[Nx] >= π) && (uphase[1] -= 2π)
        (uphase[1] - uphase[Nx] <= -π) && (uphase[1] += 2π)
end
