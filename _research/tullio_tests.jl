using Pkg
pkg"activate ."

## 
using VortexDistributions
import VortexDistributions:phase_jumps

"""
    jumps = phase_jumps(ϕ,dim)

Count phase jumps greater than `π` in phase `ϕ` along dimension `dim`.
See also: [`unwrap`](@ref), [`unwrap`](@ref), [`unwrap!`](@ref)
"""
function phase_jumps(phase,dim=1)
    @assert (dim==1 || dim==2)
    s1,s2 = size(phase)
    pdiff = zero(phase)

    if dim == 1
    for j in 1:s2
        @inbounds dϕ = phase[1,j] - phase[s1,j]
        @inbounds abs(dϕ) > π && (pdiff[1,j] += sign(dϕ))
        for i in 2:s1
            @inbounds dϕ = phase[i,j] - phase[i-1,j]
            @inbounds abs(dϕ) > π && (pdiff[i,j] += sign(dϕ))
        end

    end

    elseif dim == 2
    for i in 1:s1
        @inbounds dϕ = phase[i,1] - phase[i,s2]
        @inbounds abs(dϕ) > π && (pdiff[i,1] += sign(dϕ))
        for j in 2:s2
            @inbounds dϕ = phase[i,j] - phase[i,j-1]
            @inbounds abs(dϕ) > π && (pdiff[i,j] += sign(dϕ))
        end
    end
    end
  return pdiff
end