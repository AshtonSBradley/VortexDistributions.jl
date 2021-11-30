function get_dipoles(xi,yi,n,nplus,Lx,Ly)

  if isempty(xi) || isempty(yi)
    @info "Input position vector is empty"
    xd = []
    yd = []
    is_dipole = []
    num_found_in_dipoles = 0
    return xd,yd,is_dipole,num_found_in_dipoles
  end

  κi = [ones(nplus); -ones(n-nplus)]

  # x and y distances
  xij,yij = distances(xi,yi)

  # incorporate periodicity
  periodic_distances!(xij,yij,Lx,Ly)

  # inter-vortex distance matrix, set self distance to much larger than box size
  rij = hypot.(xij,yij)
  rij += diagm(0=>1e6*hypot(Lx,Ly)*ones(n))

  # get index of minimum separations
  _, nn_index = findmin(rij,dims=2)

  nn_mat = ones(n)*collect(1:n)' .|> Int
  nn_vec = nn_mat[nn_index][:,1]

  mutuals = (collect(1:n) .== nn_vec[nn_vec])

  # find dipoles
  is_dipole = zeros(n)
  for i in eachindex(mutuals)
    is_dipole[i] = (mutuals[i] == 1 && κi[i]*κi[nn_vec[i]] == -1.0)
  end

  dipoles_index = findall(is_dipole .> 0.0)
  dipoles_index2 = nn_vec[dipoles_index]
  num_found_in_dipoles = length(dipoles_index)

  xd = xi[dipoles_index]
  yd = yi[dipoles_index]
  xd2 = xi[dipoles_index2]
  yd2 = yi[dipoles_index2]
  xd = [xd[1:Int(end/2)]; xd2[1:Int(end/2)]]
  yd = [yd[1:Int(end/2)]; yd2[1:Int(end/2)]]

  return xd,yd,is_dipole,num_found_in_dipoles
end
