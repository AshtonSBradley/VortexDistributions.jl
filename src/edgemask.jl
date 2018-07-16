function edgemask(psi,x,y)

psi[:,1] = zeros(x)
psi[:,end] = zeros(x)
psi[1,:] = zeros(y')
psi[end,:] = zeros(y')

return psi
end
