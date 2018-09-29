#some helpers
linspace(a,b,n) = LinRange(a,b,n) |> collect

function vortexcore(r,ξ)
    return r/sqrt(r^2 + ξ^2)
end

function edgemask(psi,x,y)
    psi[:,1] = zero(x)
    psi[:,end] = zero(x)
    psi[1,:] = zero(y')
    psi[end,:] = zero(y')
return psi
end

function circmask(psi,x,y,R)
    for j in eachindex(x), k in eachindex(y)
        (x[j]^2+y[k]^2 > R^2) && (psi[j,k] = complex(0.))
    end
    return psi
end

function circmask!(phi,psi,x,y,R)
    phi .= psi
    for j in eachindex(x), k in eachindex(y)
            (x[j]^2+y[k]^2 > R^2) && (phi[j,k] = complex(0.))
    end
end

function findnz(A)
I = findall(!iszero,A)
v = A[I]
ix = [I[i][1] for i in eachindex(I)]
iy = [I[i][2] for i in eachindex(I)]
return ix,iy,v
end

function place_randomvortices(x,y,Nv)
Lx = x[end]-x[1]; Ly = x[end] - x[1]
dx = diff(x)[1]; dy = diff(y)[1]
testvort = zeros(Nv,3)
#makes sure vortices are away from edges
k = 1
while k<=Nv
a = -Lx/2 + Lx*rand()
b = -Ly/2 + Ly*rand()
σ = rand([-1,1],1)
    if (-Lx/2 + dx < a < Lx/2 - dx && -Ly/2 + dy < b < Ly/2 - dy)
        testvort[k,:] = [a b σ]
        k+=1
    end
end

testvort = sortslices(testvort,dims=1)
return testvort
end
