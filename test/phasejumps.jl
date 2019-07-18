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
