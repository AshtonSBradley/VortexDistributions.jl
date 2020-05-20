Δ(x) = x[2]-x[1]

function find_where(A)
    I = findall(!iszero,A)
    v = A[I]
    ix = [I[i][1] for i in eachindex(I)]
    iy = [I[i][2] for i in eachindex(I)]
    return ix,iy,v
end

function found_near(n)
    near = true
    for j in 1:n
        psi,vort = rand_vortexfield(1)
        vortfound = find_vortices(psi)
        vfdata = vortex_array(vortfound)
        vdata = vortex_array(vort)
        dx = Δ(psi.x)
        near *= isapprox(vdata,vfdata,rtol = dx/4)
    end
    return near
end

"""
    vortices = remove_edge_vortices(vort::Array{PointVortex,1},x,y,edge=1)

Strip artifact edgevortices arising from periodic phase differencing.
"""
function remove_edge_vortices(vort::Array{PointVortex,1},psi::Field,edge=1)
    @unpack x,y = psi; dx,dy=Δ(x),Δ(y)
    keep = []
    for j = 1:length(vort)
        xi,yi,qi = vortex_array(vort[j])
        xedge = isapprox(xi,x[1],atol=edge*dx) || isapprox(xi,x[end],atol=edge*dx)
        yedge = isapprox(yi,y[1],atol=edge*dy) || isapprox(yi,y[end],atol=edge*dy)
        not_edge = !(xedge || yedge)
        not_edge && push!(keep,j)
    end
    return vort[keep]
end


## phase jumps, unwrap

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

function phase_jumps!(pdiff,phase,dim=1)
    @assert (dim==1 || dim==2)
    s1,s2 = size(phase)

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
end

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
