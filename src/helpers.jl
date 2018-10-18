#some helpers
linspace(a,b,n) = LinRange(a,b,n) |> collect

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

function findnonzero(A)
    I = findall(!iszero,A)
    v = A[I]
    ix = [I[i][1] for i in eachindex(I)]
    iy = [I[i][2] for i in eachindex(I)]
    return ix,iy,v
end

function isinterior(a,b,x,y)
    Lx = x[end]-x[1]; Ly = y[end] - y[1]
    dx = x[2] - x[1]; dy = y[2] - y[1]
    return (-Lx/2 + 2dx < a < Lx/2 - 2dx && -Ly/2 + 2dy < b < Ly/2 - 2dy)
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

function checkvortexlocations(testvort,vortices,x,y,Nv)
#check detection to 2 x grid resolution
    dx = x[2] - x[1]; dy = y[2] - y[1]
    vortfound = 0
    for j=1:Nv
        if (isapprox(testvort[j,1],vortices[j,1],atol=2dx) && isapprox(testvort[j,2],vortices[j,2],atol=2dy))
        vortfound+=1
        end
    end
    return vortfound
end

 """
`x, DM =  chebdif(N,M)`

computes the differentiation matrices `D1, D2, ..., DM` on Chebyshev nodes.

 Input:

 `N`         Size of differentiation matrix.

 `M`         Number of derivatives required (integer).

 Note      `0 < M <= N-1`.

 Output:

 DM        `DM(1:N,1:N,ell)` contains ell-th derivative matrix, ell=1..M.

 The code implements two strategies for enhanced
 accuracy suggested by W. Don and S. Solomonoff in
 SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268 (1994).
 The two strategies are (a) the use of trigonometric
 identities to avoid the computation of differences
 x(k)-x(j) and (b) the use of the "flipping trick"
 which is necessary since sin t can be computed to high
 relative precision when t is small whereas sin (pi-t) cannot.
 Note added May 2003:  It may, in fact, be slightly better not to
 implement the strategies (a) and (b).   Please consult the following
 paper for details:   "Spectral Differencing with a Twist", by
 R. Baltensperger and M.R. Trummer, to appear in SIAM J. Sci. Comp.

 J.A.C. Weideman, S.C. Reddy 1998.  Help notes modified by
 JACW, May 2003."""

function chebdif(N, M)
    @assert 0<M<=N-1

    n1,n2 = (floor(N/2)|> Int,ceil(N/2)|> Int)  # Indices used for flipping trick.

    k = collect(0:N-1)' # Compute theta vector.
    th = k*π/(N-1)

    x = sin.(π*collect(N-1:-2:1-N)'./(2*(N-1))) # Compute Chebyshev points.

    T = repeat(th/2,N,1)'
    DX = 2*sin.(T'.+T).*sin.(T'.-T) # Trigonometric identity.
    DX = [DX[1:n1,:]; -reverse(reverse(DX[1:n2,:],dims=2),dims=1)] # Flipping trick.
    DX[diagind(DX)] .= 1 # Put 1's on the main diagonal of DX.

    C = Matrix(Toeplitz((-1).^k',(-1).^k'))  # C is the matrix with
    C[1,:] = C[1,:]*2; C[N,:] = C[N,:]*2  # entries c(k)/c(j)
    C[:,1] = C[:,1]/2; C[:,N] = C[:,N]/2

    Z = 1 ./DX  # Z contains entries 1/(x(k)-x(j))
    Z[diagind(Z)] .= 0  # with zeros on the diagonal.

    D = Matrix{Float64}(I, N, N)  # D contains diff. matrices.

    DM = zeros(N,N,M)

    for ell = 1:M
        D = ell*Z.*(C.*repeat(diag(D),1,N) - D) # Off-diagonals
        D[diagind(D)] = -sum(D',dims=1)  # Correct main diagonal of D
        DM[:,:,ell] = D  # Store current D in DM
    end

    return x, DM
    end

function getChebDMatrix(n)
    z, M = chebdif(n, 1)
    Dz = M[:,:,1]
    return z, Dz
end

function getChebD2Matrix(n)
    z, M = chebdif(n, 2)
    D2z = M[:,:,2]
    return z, D2z
end
