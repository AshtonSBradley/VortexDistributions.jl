Δ(x) = x[2]-x[1]

function find_where(A)
    I = findall(!iszero,A)
    v = A[I]
    ix = [I[i][1] for i in eachindex(I)]
    iy = [I[i][2] for i in eachindex(I)]
    return ix,iy,v
end

function found_near(n)
    print("HERE")
    near = true
    for j in 1:n
        psi,vort = rand_vortexfield(1)
        vortfound = findvortices(psi)
        vfdata = vortex_array(vortfound)
        vdata = vortex_array(vort)
        dx = Δ(psi.x)
        near *= isapprox(vdata,vfdata,rtol = dx/4)
    end
    return near
end

function found_near(n, nvorts)
    near = true
    for j in 1:n
        psi,vort = rand_vortexfield(nvorts)
        vortfound = findvortices(psi)
        vfdata = vortex_array(vortfound)
        vdata = vortex_array(vort)
        if length(vfdata[:, 1]) != nvorts
            print(length(vfdata[:, 1]))
            return false
        end
        dx = Δ(psi.x)
        sort!(vdata, dims=1)
        sort!(vfdata, dims=1)
        for i in 1:nvorts
            near *= isapprox(vdata[i, :], vfdata[i, :], rtol = dx/4)
        end
    end
    return near
end

"""
    vortices = remove_vortices_edge(vort::Array{PointVortex,1},x,y,edge=1)

Strip artifact edgevortices arising from periodic phase differencing.
"""
function remove_vortices_edge(vort::Array{PointVortex,1},psi::Field,edge=1)
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
    for j in 2:s2
        for i in 1:s1
            j==2 && (@inbounds dϕ = phase[i,1] - phase[i,s2];
            @inbounds abs(dϕ) > π && (pdiff[i,1] += sign(dϕ)))
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
    for j in 2:s2
        for i in 1:s1
            j==2 && (@inbounds dϕ = phase[i,1] - phase[i,s2];
            @inbounds abs(dϕ) > π && (pdiff[i,1] += sign(dϕ)))
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
        @inbounds (k=uphase[1,j] - uphase[s1,j];k > π) && (uphase[1,j] -= 2π*div(k,2pi,RoundNearest))
        @inbounds (k=uphase[1,j] - uphase[s1,j];k < -π) && (uphase[1,j] += 2π*div(-k,2pi,RoundNearest))
        for i in 2:s1
        @inbounds (k=uphase[i,j] - uphase[i-1,j];k > π) && (uphase[i,j] -= 2π*div(k,2pi,RoundNearest))
        @inbounds (k=uphase[i,j] - uphase[i-1,j];k < -π) && (uphase[i,j] += 2π*div(-k,2pi,RoundNearest))
        end
    end

    elseif dim == 2
        for j in 2:s2
        for i in 1:s1
            j==2 && (@inbounds (k= uphase[i,1] - uphase[i,s2]; k > π) && (uphase[i,1] -= 2π*div(k,2pi,RoundNearest));
            @inbounds (k=uphase[i,1] - uphase[i,s2];k < -π) && (uphase[i,1] += 2π*div(-k,2pi,RoundNearest)))
            @inbounds (k=uphase[i,j] - uphase[i,j-1];k > π) && (uphase[i,j] -= 2π*div(k,2pi,RoundNearest))
            @inbounds (k=uphase[i,j] - uphase[i,j-1];k < -π) && (uphase[i,j] += 2π*div(-k,2pi,RoundNearest))
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
            @inbounds (k=uphase[1,j] - uphase[s1,j];k > π) && (uphase[1,j] -= 2π*div(k,2pi,RoundNearest))
            @inbounds (k=uphase[1,j] - uphase[s1,j];k < -π) && (uphase[1,j] += 2π*div(-k,2pi,RoundNearest))
            for i in 2:s1
            @inbounds (k=uphase[i,j] - uphase[i-1,j];k > π) && (uphase[i,j] -= 2π*div(k,2pi,RoundNearest))
            @inbounds (k=uphase[i,j] - uphase[i-1,j];k < -π) && (uphase[i,j] += 2π*div(-k,2pi,RoundNearest))
            end
        end
    
        elseif dim == 2

        for j in 2:s2
            for i in 1:s1
                j==2 && (@inbounds (k= uphase[i,1] - uphase[i,s2]; k > π) && (uphase[i,1] -= 2π*div(k,2pi,RoundNearest));
                @inbounds (k=uphase[i,1] - uphase[i,s2];k < -π) && (uphase[i,1] += 2π*div(-k,2pi,RoundNearest)))
                @inbounds (k=uphase[i,j] - uphase[i,j-1];k > π) && (uphase[i,j] -= 2π*div(k,2pi,RoundNearest))
                @inbounds (k=uphase[i,j] - uphase[i,j-1];k < -π) && (uphase[i,j] += 2π*div(-k,2pi,RoundNearest))
            end
        end
    end
end

function unwrap(phase::Array{Float64,1})
    uphase = copy(phase)
    s1 = length(phase)
        @inbounds (k=uphase[1] - uphase[s1];k > π) && (uphase[1] -= 2π*div(k,2pi,RoundNearest))
        @inbounds (k=uphase[1] - uphase[s1];k < -π) && (uphase[1] += 2π*div(-k,2pi,RoundNearest))
        for i in 2:s1
        @inbounds (k=uphase[i] - uphase[i-1];k > π) && (uphase[i] -= 2π*div(k,2pi,RoundNearest))
        @inbounds (k=uphase[i] - uphase[i-1];k < -π) && (uphase[i] += 2π*div(-k,2pi,RoundNearest))
        end

        return uphase
end

function unwrap!(uphase::Array{Float64,1},phase::Array{Float64,1})
    uphase .= phase
    s1 = length(phase)
    @inbounds (k=uphase[1] - uphase[s1];k > π) && (uphase[1] -= 2π*div(k,2pi,RoundNearest))
    @inbounds (k=uphase[1] - uphase[s1];k < -π) && (uphase[1] += 2π*div(-k,2pi,RoundNearest))
    for i in 2:s1
    @inbounds (k=uphase[i] - uphase[i-1];k > π) && (uphase[i] -= 2π*div(k,2pi,RoundNearest))
    @inbounds (k=uphase[i] - uphase[i-1];k < -π) && (uphase[i] += 2π*div(-k,2pi,RoundNearest))
    end
end

function distances(x,y)
    xij = @. x - x'
    yij = @. y - y'
    return xij,yij
end

function periodic_distances!(x,y,Lx,Ly)
    make_periodic!(x,Lx)
    make_periodic!(y,Ly)
end

function make_periodic!(x,L)
    @. x[x>L/2]  -= L
    @. x[x<-L/2] += L
end

"""
    rij = sparse_distances(x,y)

Sparse form of unique lower triangular distance matrix for vectors `x`, `y` specifying points on the plane.
"""
function sparse_distances(x,y)
    xij = @. x - x'
    yij = @. y - y'
    rij = hypot.(xij,yij) |> tril |> sparse
end

"""
    tree = spanning_tree(x,y)

Minimal spanning tree for points specified by vectors `x`, `y`.
"""
function spanning_tree(x,y)
    rij = sparse_distances(x,y)
    # graph representation
    g = SimpleWeightedGraph(length(rij))
    for edge in eachindex(rij)
        n1,n2 = edge[1],edge[2]
        add_edge!(g,n1,n2,rij[n1,n2])
    end
    # get the spanning tree
    tree = prim_mst(g)
end

"""
    vortex_marker(c)

Define vortex marker using sign.
"""
function vortex_marker(c)
    if c==1.0
        return (4,0.6,:green, stroke(1.3, 1.0, :white))
    elseif c==-1.0
        return (4,0.6,:blue, stroke(1.3, 1.0, :white))
    else
    end
end

# """
#     plot_vortices!(p,vdata)

# Plot vortices defined in `vdata` on figure `p`.
# """
# function plot_vortices!(p,vdata)
#     for j in 1:size(vdata)[1]
#         vx,vy = vdata[j,1],vdata[j,2]
#         scatter!(p1,[vx],[vy],marker=vortex_marker(vdata[j,3]),label=false)
#     end
#     return p1
# end

# """
#     plot_cluster_vortices!(p,cluster::Cluster)

# Plot cluster in figure `p` using `vortex_marker`.
# """
# function plot_cluster_vortices!(p,cluster::Cluster)
#     vdata = vortex_array(cluster.vortices)
#     for j in 1:size(vdata)[1]
#         vx,vy = vdata[j,1],vdata[j,2]
#         scatter!(p,[vx],[vy],marker=vortex_marker(vdata[j,3]),label=false)
#     end
#     return p
# end

# """
#     plot_cluster_tree!(p,cluster::Cluster)

# Plot minimal spanning tree on figure `p` for `cluster`.
# """
# function plot_cluster_tree!(p,cluster::Cluster)
#     vortices = cluster.vortices
#     tree = cluster.tree
#     va = vortex_array(vortices)
#     xi = @view va[:,1]
#     yi = @view va[:,2]
#     ci = @view va[:,3]

#     for i in 1:length(tree)
#         edge = tree[i]
#         src,dst = edge.src,edge.dst
#         x1,y1 = xi[src],yi[src]
#         x2,y2 = xi[dst],yi[dst]
#         csign = vortices[1].qv
#         col = csign == 1.0 ? :green : :blue
#         plot!(p,[x1,x2],[y1,y2],label=false,c=col,alpha=0.4)
#     end
#     return p
# end

# """
#     plot_custer!(p,c)

# Plot the vortices and minimal spanning tree.
# """
# function plot_cluster!(p,c)
#     plot_cluster_vortices!(p,c)
#     plot_cluster_tree!(p,c)
#     return p
# end
