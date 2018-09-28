#some helpers
linspace(a,b,n) = LinRange(a,b,n) |> collect

function vortexcore(r,ξ)
    return r/sqrt(r^2 + ξ^2)
end

function edgemask(psi,x,y)
    psi[:,1] = zeros(x)
    psi[:,end] = zeros(x)
    psi[1,:] = zeros(y')
    psi[end,:] = zeros(y')
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
