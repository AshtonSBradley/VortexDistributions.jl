"""
jumps = countphasejumps(phase,dim)

Count jumps greater than π in `phase` along dimension `dim`
"""
function countphasejumps(phase,dim=1)
    @assert (dim==1 || dim==2)
    Nx,Ny = size(phase)
    pdiff = zeros(phase)

    if dim == 1
    for j in 1:Ny
        for i in 2:Nx
            Δ = phase[i,j] - phase[i-1,j]
            abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end
            Δ = phase[1,j] - phase[Nx,j]
            abs(Δ) > π && (pdiff[1,j] += sign(Δ))
    end

    elseif dim == 2
    for j in 2:Ny
        for i in 1:Nx
            Δ = phase[i,j] - phase[i,j-1]
            abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end
    end
        for i in 1:Nx
            Δ = phase[i,1] - phase[i,Ny]
            abs(Δ) > π && (pdiff[i,1] += sign(Δ))
        end
    end

  return pdiff
end

function countphasejumps!(pdiff,phase,dim=1)
    @assert (dim==1 || dim==2)
    Nx,Ny = size(phase)

    if dim == 1
    for j in 1:Ny
        for i in 2:Nx
            Δ = phase[i,j] - phase[i-1,j]
            abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end
            Δ = phase[1,j] - phase[Nx,j]
            abs(Δ) > π && (pdiff[1,j] += sign(Δ))
    end

    elseif dim == 2
    for j in 2:Ny
        for i in 1:Nx
            Δ = phase[i,j] - phase[i,j-1]
            abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end

    end
        for i in 1:Nx
            Δ = phase[i,1] - phase[i,Ny]
            abs(Δ) > π && (pdiff[i,1] += sign(Δ))
        end
    end
end
