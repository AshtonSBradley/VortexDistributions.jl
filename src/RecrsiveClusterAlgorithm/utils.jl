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

"""
    plot_vortices!(p,vdata)

Plot vortices defined in `vdata` on figure `p`.
"""
function plot_vortices!(p,vdata)
    for j in 1:size(vdata)[1]
        vx,vy = vdata[j,1],vdata[j,2]
        scatter!(p1,[vx],[vy],marker=vortex_marker(vdata[j,3]),label=false)
    end
    return p1
end

"""
    plot_cluster_vortices!(p,cluster::Cluster)

Plot cluster in figure `p` using `vortex_marker`.
"""
function plot_cluster_vortices!(p,cluster::Cluster)
    vdata = vortex_array(cluster.vortices)
    for j in 1:size(vdata)[1]
        vx,vy = vdata[j,1],vdata[j,2]
        scatter!(p,[vx],[vy],marker=vortex_marker(vdata[j,3]),label=false)
    end
    return p
end

"""
    plot_cluster_tree!(p,cluster::Cluster)

Plot minimal spanning tree on figure `p` for `cluster`.
"""
function plot_cluster_tree!(p,cluster::Cluster)
    vortices = cluster.vortices
    tree = cluster.tree
    va = vortex_array(vortices)
    xi = @view va[:,1]
    yi = @view va[:,2]
    ci = @view va[:,3]

    for i in 1:length(tree)
        edge = tree[i]
        src,dst = edge.src,edge.dst
        x1,y1 = xi[src],yi[src]
        x2,y2 = xi[dst],yi[dst]
        csign = vortices[1].qv
        col = csign == 1.0 ? :green : :blue
        plot!(p,[x1,x2],[y1,y2],label=false,c=col,alpha=0.4)
    end
    return p
end

"""
    plot_custer!(p,c)

Plot the vortices and minimal spanning tree.
"""
function plot_cluster!(p,c)
    plot_cluster_vortices!(p,c)
    plot_cluster_tree!(p,c)
    return p
end
