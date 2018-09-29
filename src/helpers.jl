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

function isinterior(a,b,x,y)
    Lx = x[end]-x[1]; Ly = y[end] - y[1]
    dx = x[2] - x[1]; dy = y[2] - y[1]
    return (-Lx/2 + dx < a < Lx/2 - dx && -Ly/2 + dy < b < Ly/2 - dy)
end

function randomvortices(x,y,Nv)
    Lx = x[end]-x[1]; Ly = y[end] - y[1]
    dx = x[2] - x[1]; dy = y[2] - y[1]
    testvort = zeros(Nv,3)

    k = 1
    while k<=Nv
        a = -Lx/2 + Lx*rand()
        b = -Ly/2 + Ly*rand()
        σ = rand([-1 1])

        #make sure vortices are away from edges
        if isinterior(a,b,x,y)
            testvort[k,:] = [a b σ]
            k+=1
        end
    end
    return sortslices(testvort,dims=1)
end
