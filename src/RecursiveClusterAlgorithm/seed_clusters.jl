function seed_clusters(xi,yi,n,nplus,Lx,Ly)

    κi = [ones(nplus); -ones(n-nplus)]

    # x and y distances
    xij = xi .- xi'
    yij = yi .- yi'

    # periodicity
   periodic_distances!(xij,yij,Lx,Ly)

    # distances: set self distance to much larger than box size
    rij = @. hypot(xij,yij)
    rij += diagm(0=>1e6*hypot(Lx,Ly)*ones(n))

    # get index of minimum separations
    _, nn_index = findmin(rij,dims=2)

    nn_mat = ones(n)*collect(1:n)' .|> Int
    nn_vec = nn_mat[nn_index][:,1]
    mutuals = (collect(1:n) .== nn_vec[nn_vec])

    is_seed = zeros(n)
    for i in eachindex(mutuals)
        is_seed[i] = (mutuals[i] == 1 && κi[i]*κi[nn_vec[i]] == 1.0)
    end

    seeds_index = findall(is_seed .> 0.0)
    seeds_index2= nn_vec[seeds_index]
    num_found_in_seeds = length(seeds_index)

    # seed coords
    xs = xi[seeds_index]
    ys = yi[seeds_index]
    xs2 = xi[seeds_index2]
    ys2 = yi[seeds_index2]

    xs = [xs; xs2]
    ys = [ys; ys2]
    seeds_index = [seeds_index; seeds_index2]

return xs,ys,is_seed,seeds_index,num_found_in_seeds
end
