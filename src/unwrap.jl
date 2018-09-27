"""
unwrapped = unwrap(phase,dim=1)

Unwraps 2d array `phase` along dimension `dim`.

unwrap!(unwrapped,phase,dim) writes in-place to unwrapped.
"""
function unwrap(phase::Array{Float64,2},dim=1)
    @assert (dim==1 || dim==2)
    Nx,Ny = size(phase)
    uphase = copy(phase)

    if dim == 1
    for j in 1:Ny
        for i in 2:Nx
        (uphase[i,j] - uphase[i-1,j] >= π) && (uphase[i,j] -= 2π)
        (uphase[i,j] - uphase[i-1,j] <= -π) && (uphase[i,j] += 2π)
        end
        (uphase[1,j] - uphase[Nx,j] >= π) && (uphase[1,j] -= 2π)
        (uphase[1,j] - uphase[Nx,j] <= -π) && (uphase[1,j] += 2π)
    end

    elseif dim == 2
    for j in 2:Ny
        for i in 1:Nx
        (uphase[i,j] - uphase[i,j-1] >= π) && (uphase[i,j] -= 2π)
        (uphase[i,j] - uphase[i,j-1] <= -π) && (uphase[i,j] += 2π)
        end

    end
        for i in 1:Nx
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
    for j in 1:Ny
        for i in 2:Nx
        (uphase[i,j] - uphase[i-1,j] >= π) && (uphase[i,j] -= 2π)
        (uphase[i,j] - uphase[i-1,j] <= -π) && (uphase[i,j] += 2π)
        end
        (uphase[1,j] - uphase[Nx,j] >= π) && (uphase[1,j] -= 2π)
        (uphase[1,j] - uphase[Nx,j] <= -π) && (uphase[1,j] += 2π)
    end

    elseif dim == 2
    for j in 2:Ny
        for i in 1:Nx
        (uphase[i,j] - uphase[i,j-1] >= π) && (uphase[i,j] -= 2π)
        (uphase[i,j] - uphase[i,j-1] <= -π) && (uphase[i,j] += 2π)
        end

    end
        for i in 1:Nx
        (uphase[i,1] - uphase[i,Ny] >= π) && (uphase[i,1] -= 2π)
        (uphase[i,1] - uphase[i,Ny] <= -π) && (uphase[i,1] += 2π)
        end
    end
end

function unwrap(phase::Array{Float64,1})
    uphase = copy(phase)
    Nx = length(phase)
        for i in 2:Nx
        (uphase[i] - uphase[i-1] >= π) && (uphase[i] -= 2π)
        (uphase[i] - uphase[i-1] <= -π) && (uphase[i] += 2π)
        end
        (uphase[1] - uphase[Nx] >= π) && (uphase[1] -= 2π)
        (uphase[1] - uphase[Nx] <= -π) && (uphase[1] += 2π)
        return uphase
end

function unwrap!(uphase::Array{Float64,1},phase::Array{Float64,1})
    Nx = length(phase)
        for i in 2:Nx
        (uphase[i] - uphase[i-1] >= π) && (uphase[i] -= 2π)
        (uphase[i] - uphase[i-1] <= -π) && (uphase[i] += 2π)
        end
        (uphase[1] - uphase[Nx] >= π) && (uphase[1] -= 2π)
        (uphase[1] - uphase[Nx] <= -π) && (uphase[1] += 2π)
end
