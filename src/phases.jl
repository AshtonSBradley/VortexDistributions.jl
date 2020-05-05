# phase jumps, unwrap
"""
    ϕu = unwrap(ϕ,dim=1)

Unwraps 2d array `ϕ` along dimension `dim`, acting periodically to give back array of same size as `ϕ`.

See also: [`unwrap!`](@ref)
"""
function unwrap(phase::Array{Float64,2},dim=1)
    @assert (dim==1 || dim==2)
    s1,s2 = size(phase)
    uphase = copy(phase)

    if dim == 1
    for j in 1:s2
        @inbounds uphase[1,j] - uphase[s1,j] >= π && (uphase[1,j] -= 2π)
        @inbounds uphase[1,j] - uphase[s1,j] <= -π && (uphase[1,j] += 2π)
        for i in 2:s1
        @inbounds uphase[i,j] - uphase[i-1,j] >= π && (uphase[i,j] -= 2π)
        @inbounds uphase[i,j] - uphase[i-1,j] <= -π && (uphase[i,j] += 2π)
        end
    end

    elseif dim == 2
    for i in 1:s1
            @inbounds uphase[i,1] - uphase[i,s2] >= π && (uphase[i,1] -= 2π)
            @inbounds uphase[i,1] - uphase[i,s2] <= -π && (uphase[i,1] += 2π)
        for j in 2:s2
            @inbounds uphase[i,j] - uphase[i,j-1] >= π && (uphase[i,j] -= 2π)
            @inbounds uphase[i,j] - uphase[i,j-1] <= -π && (uphase[i,j] += 2π)
        end
    end
end

    return uphase
end

"""
    unwrap!(ϕu,ϕ,dim)

Unwraps 2d phase array `ϕ` along dimension `dim`,
acting periodically and writing the unwrapped array `ϕ`
in place.

See also: [`unwrap`](@ref)
"""
function unwrap!(uphase::Array{Float64,2},phase::Array{Float64,2},dim=1)
    @assert (dim==1 || dim==2)
    s1,s2 = size(phase)
    uphase .= phase

    if dim == 1
    for j in 1:s2
        @inbounds uphase[1,j] - uphase[s1,j] >= π && (uphase[1,j] -= 2π)
        @inbounds uphase[1,j] - uphase[s1,j] <= -π && (uphase[1,j] += 2π)
        for i in 2:s1
        @inbounds uphase[i,j] - uphase[i-1,j] >= π && (uphase[i,j] -= 2π)
        @inbounds uphase[i,j] - uphase[i-1,j] <= -π && (uphase[i,j] += 2π)
        end
    end

    elseif dim == 2
    for i in 1:s1
            @inbounds uphase[i,1] - uphase[i,s2] >= π && (uphase[i,1] -= 2π)
            @inbounds uphase[i,1] - uphase[i,s2] <= -π && (uphase[i,1] += 2π)
        for j in 2:s2
            @inbounds uphase[i,j] - uphase[i,j-1] >= π && (uphase[i,j] -= 2π)
            @inbounds uphase[i,j] - uphase[i,j-1] <= -π && (uphase[i,j] += 2π)
        end
    end
end
end

function unwrap(phase::Array{Float64,1})
    uphase = copy(phase)
    s1 = length(phase)
        @inbounds uphase[1] - uphase[s1] >= π && (uphase[1] -= 2π)
        @inbounds uphase[1] - uphase[s1] <= -π && (uphase[1] += 2π)
        for i in 2:s1
        @inbounds uphase[i] - uphase[i-1] >= π && (uphase[i] -= 2π)
        @inbounds uphase[i] - uphase[i-1] <= -π && (uphase[i] += 2π)
        end

        return uphase
end

function unwrap!(uphase::Array{Float64,1},phase::Array{Float64,1})
    s1 = length(phase)
        @inbounds uphase[1] - uphase[s1] >= π && (uphase[1] -= 2π)
        @inbounds uphase[1] - uphase[s1] <= -π && (uphase[1] += 2π)
        for i in 2:s1
        @inbounds uphase[i] - uphase[i-1] >= π && (uphase[i] -= 2π)
        @inbounds uphase[i] - uphase[i-1] <= -π && (uphase[i] += 2π)
        end
end

"""
    jumps = phasejumps(ϕ,dim)

Count phase jumps greater than `π` in phase `ϕ` along dimension `dim`.

See also: [`unwrap`](@ref), [`unwrap`](@ref), [`unwrap!`](@ref)
"""
function phasejumps(phase,dim=1)
    @assert (dim==1 || dim==2)
    s1,s2 = size(phase)
    pdiff = zero(phase)

    if dim == 1
    for j in 1:s2
        @inbounds Δ = phase[1,j] - phase[s1,j]
        @inbounds abs(Δ) > π && (pdiff[1,j] += sign(Δ))
        for i in 2:s1
            @inbounds Δ = phase[i,j] - phase[i-1,j]
            @inbounds abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end

    end

    elseif dim == 2
    for i in 1:s1
        @inbounds Δ = phase[i,1] - phase[i,s2]
        @inbounds abs(Δ) > π && (pdiff[i,1] += sign(Δ))
        for j in 2:s2
            @inbounds Δ = phase[i,j] - phase[i,j-1]
            @inbounds abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end
    end
    end
  return pdiff
end

function phasejumps!(pdiff,phase,dim=1)
    @assert (dim==1 || dim==2)
    s1,s2 = size(phase)

    if dim == 1
    for j in 1:s2
        @inbounds Δ = phase[1,j] - phase[s1,j]
        @inbounds abs(Δ) > π && (pdiff[1,j] += sign(Δ))
        for i in 2:s1
            @inbounds Δ = phase[i,j] - phase[i-1,j]
            @inbounds abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end

    end

    elseif dim == 2
    for i in 1:s1
        @inbounds Δ = phase[i,1] - phase[i,s2]
        @inbounds abs(Δ) > π && (pdiff[i,1] += sign(Δ))
        for j in 2:s2
            @inbounds Δ = phase[i,j] - phase[i,j-1]
            @inbounds abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end
    end
    end
end
