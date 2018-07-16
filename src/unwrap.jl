"""

unwrap(phase,dim=1)

Given an input phase `phase` and dimension to unwrap `dim`, returns an unwrapped array.

unwrap!(unwrapped,phase,dim) writes in-place to unwrapped. 
"""
function unwrap(phase,dim=1)
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

function unwrap!(uphase,phase,dim=1)
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

#= old
function unwrap(phase, inplace=false)

  unwrapped = inplace ? phase : copy(phase)
  for i in 2:length(phase)
    unwrapped[i] - unwrapped[i-1] >= π && (unwrapped[i] -= 2π)
    unwrapped[i] - unwrapped[i-1] <= -π && (unwrapped[i] += 2π)
  end
  return unwrapped
end

unwrap!(phase) = unwrap(phase, true)
=#
