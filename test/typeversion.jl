using Parameters, Interpolations

abstract type Topology end

struct Torus <: Topology
    x::Vector{Float64}
    y::Vector{Float64}
    ψ::Array{Complex{Float64},2}
end

struct Sphere <: Topology
    x::Vector{Float64}
    y::Vector{Float64}
    ψ::Array{Complex{Float64},2}
end

"""
    ix,iy,v = findwhere(A)
Finds indices `ix`, `iy`, and values where `A` returns true.
"""
function findwhere(A)
    I = findall(!iszero,A)
    v = A[I]
    ix = [I[i][1] for i in eachindex(I)]
    iy = [I[i][2] for i in eachindex(I)]
    return ix,iy,v
end

"""
    jumps = phasejumps(phase,dim)

Count phase jumps greater than π in `phase` along dimension `dim`
"""
function phasejumps(phase,dim=1)
    @assert (dim==1 || dim==2)
    Nx,Ny = size(phase)
    pdiff = zero(phase)

    if dim == 1
    @inbounds for j in 1:Ny
        @inbounds for i in 2:Nx
            Δ = phase[i,j] - phase[i-1,j]
            abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end
            Δ = phase[1,j] - phase[Nx,j]
            abs(Δ) > π && (pdiff[1,j] += sign(Δ))
    end

    elseif dim == 2
    @inbounds for j in 2:Ny
        @inbounds for i in 1:Nx
            Δ = phase[i,j] - phase[i,j-1]
            abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end
    end
        @inbounds for i in 1:Nx
            Δ = phase[i,1] - phase[i,Ny]
            abs(Δ) > π && (pdiff[i,1] += sign(Δ))
        end
    end

  return pdiff
end

function phasejumps!(pdiff,phase,dim=1)
    @assert (dim==1 || dim==2)
    Nx,Ny = size(phase)

    if dim == 1
    @inbounds for j in 1:Ny
        @inbounds for i in 2:Nx
            Δ = phase[i,j] - phase[i-1,j]
            abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end
            Δ = phase[1,j] - phase[Nx,j]
            abs(Δ) > π && (pdiff[1,j] += sign(Δ))
    end

    elseif dim == 2
    @inbounds for j in 2:Ny
        @inbounds for i in 1:Nx
            Δ = phase[i,j] - phase[i,j-1]
            abs(Δ) > π && (pdiff[i,j] += sign(Δ))
        end

    end
        @inbounds for i in 1:Nx
            Δ = phase[i,1] - phase[i,Ny]
            abs(Δ) > π && (pdiff[i,1] += sign(Δ))
        end
    end
end

"""
    `nt,np,nn,vortices = findvortices(psi<:Topology;shift=true)`
Returns `nt` total vortices, consisting of `np` positive, `nn` negative, and housed in the array `vortices`.
"""
function findvorticesgrid(psi::Torus;shift=true)
    @unpack x,y,ψ = psi
    phase = angle.(ψ)
    diffx = phasejumps(phase,1)
    diffy = phasejumps(phase,2)

    # use in-place memory recycling
    circshift!(phase,diffx,(0,1))
    diffx .-= phase
    diffx .-= diffy
    circshift!(phase,diffy,(1,0))
    diffx .+= phase

    ixp,iyp,vp = findwhere(diffx .> 0.0)
    xp = x[ixp]; yp = y[iyp]; np = length(vp)

    ixn,iyn,vn = findwhere(diffx .< 0.0)
    xn = x[ixn]; yn = y[iyn]; nn = length(vn)

    nt = np + nn

    # TODO this may need to be revisited (related to bug?)
    if shift
        dx = x[2]-x[1]; dy = y[2]-y[1]
        xp .+= -dx/2; yp .+= -dy/2; xn .+= -dx/2; yn .+= -dy/2
    end

    vortices = [xn yn -vn; xp yp vp]
    vortices = sortslices(vortices,dims=1)

    return nt,np,nn,vortices
end

function findvorticesgrid(psi::Sphere;shift=true)
    @unpack x,y,ψ = psi
    phase = angle.(ψ)
    diffx = phasejumps(phase,1)
    diffy = phasejumps(phase,2)

    #use in-place memory recycling
    circshift!(phase,diffx,(0,1))
    diffx .-= phase
    diffx .-= diffy
    circshift!(phase,diffy,(1,0))
    diffx .+= phase

    ixp,iyp,vp = findwhere(diffx.>0.)
    xp = x[ixp]; yp = y[iyp]; np = length(vp)

    ixn,iyn,vn = findwhere(diffx.<0.)
    xn = x[ixn]; yn = y[iyn]; nn = length(vn)

    nt = np + nn

    # TODO this may need to be revisited (related to bug?)
    if shift
        dx = x[2]-x[1]; dy = y[2] - y[1]
        xp .+= -dx/2; yp .+= -dy/2; xn .+= -dx/2; yn .+= -dy/2
    end

    #find polar winding
    vortices = [xn yn -vn; xp yp vp]

    # to do: optimize this step at start of method to avoid extra angle call:
    windvals = phasejumps(angle.(ψ),2)

    #test north pole for winding. - sign means polar vortex co-rotating with w > 0
    w1 = sum(windvals[1,:])
    (sign(w1) != 0) && (vortices = [vortices; 0.0 0.0 -w1])

    #test south pole for winding. + sign means co-rotating with w > 0
    w2 = sum(windvals[end,:])
    (sign(w2) != 0) && (vortices = [vortices; pi 0.0 w2])

    vortices = sortslices(vortices,dims=1)
    return nt,np,nn,vortices
end


function findvorticesinterp(psi)
    @unpack ψ,x,y = psi
    nt,np,nn,vortices = findvorticesgrid(psi;shift=false)
    nt,np,nn,vortices = removeedgevortices(vortices,x,y)
    #TODO: allow for interp with periodic data (here edges are stripped to avoid)
    for j in 1:nt
        vortex = vortices[j,:]
        vortz,psiz,xz,yz = corezoom(vortex,ψ,x,y)
        vortz,psiz,xz,yz = corezoom(vortz,psiz,xz,yz)
        vortz,psiz,xz,yz = corezoom(vortz,psiz,xz,yz)
        vortices[j,1:2] = vortz[1:2]
    end
    return nt,np,nn,vortices
end

"""
    nt,np,nn,vortices = findvortices(ψ,x,y)

Locates vortices as 2π phase windings around plaquettes on a cartesian spatial grid. Uses an optimized plaquette method followed by recursive interpolation.

Requires a 2D wavefunction ψ(x,y) on a cartesian grid specified by vectors x, y.

`nt` - total number of vortices

`np` - number of positive vortices

`nn` - number of negative vortices

`vortices` - array of vortex coordinates `xv,yv` and circulations `cv`. Each row is of the form `[xv, yv, cv]`, and the array is sorted into lexical order according to the `xv` coordinates
"""
function findvortices(psi::T;interp=true) where T<:Topology
    if interp
        return nt,np,nn,vortices = findvorticesinterp(psi)
    else
        return nt,np,nn,vortices = findvorticesgrid(psi)
    end
end

function findvortexmask(psi::T,R) where T<: Topology
    @unpack x,y,ψ = psi
    ψ = circmask(ψ,x,y,1.1*R)
    ϕ = T(x,y,ψ)
    nt,np,nn,vortices = findvortices(ϕ)

# TODO more tests
# remove vortices found outside mask boundary
for (i,xv) in enumerate(vortices[:,1]), yv in vortices[i,2]
    (norm([xv,yv]) > R) && (vortices[i,:] = [0. 0. 0.])
end

vortfound = Array{Float64, 2}[]
for (i,sv) in enumerate(vortices[:,3])
    (sv!=0) && push!(vortfound,vortices[i,:]')
end
    # currently will fail for more that one vortex (?)
return vortfound[1]
end

"""
`vortices = remove_edgevortices(vortices,x,y,edge=1)`

Strips edgevortices due to periodic phase differencing."""
function removeedgevortices(vortices,x,y,edge=1)
    dx=x[2]-x[1]; dy=y[2]-y[1]
    Nfound,_ = size(vortices)
    keep = []
    for j = 1:Nfound
        xi = vortices[j,1]; yi = vortices[j,2];
        xedge = isapprox(xi,x[1],atol=edge*dx) || isapprox(xi,x[end],atol=edge*dx)
        yedge = isapprox(yi,y[1],atol=edge*dy) || isapprox(yi,y[end],atol=edge*dy)
        not_edge = !(xedge || yedge)
        not_edge && push!(keep,j)
    end
    vortices = vortices[keep,:]
    nt = size(vortices)[1]
    np = sum(vortices[:,3].>0)
    nn = nt - np
    return nt,np,nn,vortices
end

"""
    vortz,psiz,xz,yz = corezoom(vortices,psi,x,y,winhalf=2,Nz=30)

Uses local interpolation to resolve core location to 4-5 figures.
"""
function corezoom(vortices,psi,x,y,winhalf=2,Nz=30)
    xv,yv = vortices[1:2]
    dx=x[2]-x[1];dy=y[2]-y[1]
    ixv = isapprox.(x,xv,atol=dx) |> findlast
    iyv = isapprox.(y,yv,atol=dy) |> findlast
    ixwin = (ixv-winhalf):(ixv+winhalf-1)
    iywin = (iyv-winhalf):(iyv+winhalf-1)
    xw = x[ixwin]; yw = y[iywin]; psiw = psi[ixwin,iywin]
    xz = LinRange(xw[1],xw[end],Nz)
    yz = LinRange(yw[1],yw[end],Nz)
    knots = (xw,yw)
    itp = interpolate(knots, psiw, Gridded(Linear()))
    psiz = itp(xz,yz)
    ψ = Torus(xz |> Vector,yz |> Vector,psiz)
    nt,np,nn,vortz = findvorticesgrid(ψ;shift=false)
    nt,np,nn,vortz = removeedgevortices(vortz,xz |> Vector,yz |> Vector)
    return vortz,psiz,xz,yz
end


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


function makevortex(ψ,vortex,x,y,ξ=1.0)
    @assert typeof(x)==Array{Float64,1}
    @assert typeof(y)==Array{Float64,1}
    @assert typeof(ψ)==Array{Complex{Float64},2}
    x0, y0, σ0 = vortex
    R(x,y) = sqrt(x^2+y^2)
    return  @. ψ*vortexcore(R(x.-x0,y'.-y0),ξ)*exp(im*σ0*atan(y'.-y0,x.-x0))
end

function makevortex!(ψ,vortex,x,y,ξ=1.0)
    @assert typeof(x)==Array{Float64,1}
    @assert typeof(y)==Array{Float64,1}
    @assert typeof(ψ)==Array{Complex{Float64},2}
    x0, y0, σ0 = vortex
    R(x,y) = sqrt(x^2+y^2)
    ψ .= @. ψ*vortexcore(R(x.-x0,y'.-y0),ξ)*exp(im*σ0*atan(y'.-y0,x.-x0))
end

function makeallvortices!(ψ,vortices,x,y,ξ=1.0)
    for j = 1:size(vortices)[1]
        makevortex!(ψ,vortices[j,:],x,y,ξ)
    end
end


"""
`unwrapped = unwrap(phase,dim=1)`

Unwraps 2d array `phase` along dimension `dim`, acting periodically to give back array of same size as `phase`.

`unwrap!(unwrapped,phase,dim)` writes in-place to `unwrapped`.
"""
function unwrap(phase::Array{Float64,2},dim=1)
    @assert (dim==1 || dim==2)
    Nx,Ny = size(phase)
    uphase = copy(phase)

    if dim == 1
    @inbounds for j in 1:Ny
        @inbounds for i in 2:Nx
        (uphase[i,j] - uphase[i-1,j] >= π) && (uphase[i,j] -= 2π)
        (uphase[i,j] - uphase[i-1,j] <= -π) && (uphase[i,j] += 2π)
        end
        (uphase[1,j] - uphase[Nx,j] >= π) && (uphase[1,j] -= 2π)
        (uphase[1,j] - uphase[Nx,j] <= -π) && (uphase[1,j] += 2π)
    end

    elseif dim == 2
    @inbounds for j in 2:Ny
        @inbounds for i in 1:Nx
        (uphase[i,j] - uphase[i,j-1] >= π) && (uphase[i,j] -= 2π)
        (uphase[i,j] - uphase[i,j-1] <= -π) && (uphase[i,j] += 2π)
        end

    end
        @inbounds for i in 1:Nx
        (uphase[i,1] - uphase[i,Ny] >= π) && (uphase[i,1] -= 2π)
        (uphase[i,1] - uphase[i,Ny] <= -π) && (uphase[i,1] += 2π)
        end
    end

  return uphase
end

function unwrap!(uphase::Array{Float64,2},phase::Array{Float64,2},dim=1)
    @assert (dim==1 || dim==2)
    Nx,Ny = size(phase)
    uphase .= phase

    if dim == 1
    @inbounds for j in 1:Ny
        @inbounds for i in 2:Nx
        (uphase[i,j] - uphase[i-1,j] >= π) && (uphase[i,j] -= 2π)
        (uphase[i,j] - uphase[i-1,j] <= -π) && (uphase[i,j] += 2π)
        end
        (uphase[1,j] - uphase[Nx,j] >= π) && (uphase[1,j] -= 2π)
        (uphase[1,j] - uphase[Nx,j] <= -π) && (uphase[1,j] += 2π)
    end

    elseif dim == 2
    @inbounds for j in 2:Ny
        @inbounds for i in 1:Nx
        (uphase[i,j] - uphase[i,j-1] >= π) && (uphase[i,j] -= 2π)
        (uphase[i,j] - uphase[i,j-1] <= -π) && (uphase[i,j] += 2π)
        end

    end
        @inbounds for i in 1:Nx
        (uphase[i,1] - uphase[i,Ny] >= π) && (uphase[i,1] -= 2π)
        (uphase[i,1] - uphase[i,Ny] <= -π) && (uphase[i,1] += 2π)
        end
    end
end

function unwrap(phase::Array{Float64,1})
    uphase = copy(phase)
    Nx = length(phase)
        @inbounds for i in 2:Nx
        (uphase[i] - uphase[i-1] >= π) && (uphase[i] -= 2π)
        (uphase[i] - uphase[i-1] <= -π) && (uphase[i] += 2π)
        end
        (uphase[1] - uphase[Nx] >= π) && (uphase[1] -= 2π)
        (uphase[1] - uphase[Nx] <= -π) && (uphase[1] += 2π)
        return uphase
end

function unwrap!(uphase::Array{Float64,1},phase::Array{Float64,1})
    Nx = length(phase)
        @inbounds for i in 2:Nx
        (uphase[i] - uphase[i-1] >= π) && (uphase[i] -= 2π)
        (uphase[i] - uphase[i-1] <= -π) && (uphase[i] += 2π)
        end
        (uphase[1] - uphase[Nx] >= π) && (uphase[1] -= 2π)
        (uphase[1] - uphase[Nx] <= -π) && (uphase[1] += 2π)
end


# vortex core construction

function vortexcore(r,ξ,ansatz=true,charge=1)
    ansatz ? (return r/sqrt(r^2 + ξ^2)) : (return core_chargen(r/ξ,charge))
end

"""
Make and evaluate the vortex core interpolation for charge `n`
"""
function core_chargen(x,n,ξ=1)
    y,ψ,res = gpecore(n)
    ψi = interpolate(tuple(y[1:end-1]), ψ[1:end-1], Gridded(Linear()))
    return ψi(x/ξ)
end

"""
Make the vortex core interpolation for charge `n`
"""
function make_fastcore(n)
    y,ψ,res = gpecore(n)
    ψi = interpolate(tuple(y[1:end-1]), ψ[1:end-1], Gridded(Linear()))
    return ψi
end

"""
Evaluate fast core interpolation
"""
function vortexcore(r,ψi::Interpolations.GriddedInterpolation,ξ=1)
    return ψi(r/ξ)
end

"""
Evaluate the gpe for vortex core solution. Slow: for testing and initializing.
"""
function gpecore(K,L=2,N=100,R = K)
    #currently r does  nothing!
    #N = 100
    #L = 2
    #Κ = 1
    #R = 2 # Stretching coordinate.  R ~ kappa seems to be a good ballpark

# Convert Chebyshev grid to [0,L] from [1 -1];
z,Dz = getChebDMatrix(N)
blank,D2z = getChebD2Matrix(N)
z = vec(reverse(z,dims=2))
z = (z)*L./2
Dz = -2*Dz./L
D2z = 4*D2z./L.^2

# y is the physical coordinate, range [0 infinity] (i.e y = r)
y =  R*(1 .+z)./(1 .-z)

#initial guess based on ansatz
Λ = 0.8248
ψ = @. y/sqrt(y^2 + Λ^-2)
ψ[1] = 0
ψ[end] = 1
ψ0 = vec(ψ)

Q = vec(z.^2 .- 2*z .+ 1)
Qmat = repeat(Q,1,N)
Zmat = repeat(2*(z .-1),1,N)
Ymat = repeat(1 ./y,1,N);

# Second Derivative
residuals = -0.5*( (Q.*(D2z*ψ) + 2*(z .-1).*(Dz*ψ) ).*Q/(4*R^2)
        + (Q/(2*R) ./y).*(Dz*ψ) )+ 0.5*K^2 *ψ ./y.^2 + ψ.^3 .- ψ
residuals[1] = 0; residuals[end] =0

while sum(abs.(residuals).^2) > 1e-12
    # Second Derivative
    residuals = -0.5*( (Q.*(D2z*ψ) + 2*(z .-1).*(Dz*ψ) ).*Q/(4*R^2)
            + (Q/(2*R) ./y).*(Dz*ψ) )+ 0.5*K^2 *ψ ./y.^2 + ψ.^3 .- ψ
    residuals[1] = 0; residuals[end] =0

    #println(sum(abs.(residuals).^2))


    Jacobi = -0.5*( (Qmat.*(D2z) + Zmat.*(Dz) ).*(Qmat/(4*R^2))
    + (Qmat/(2*R).*Ymat).*(Dz) ) + diagm(0 => 0.5*K^2 ./y.^2)+ diagm(0 => 3*ψ.^2 .- 1)


    Jacobi[1,:] = [1 zeros(1,N-1)]
    Jacobi[N,:] = [zeros(1,N-1) 1]


    Δ = vec(-4/7*(Jacobi\residuals))

    ψ = ψ + Δ
    ψ[1] = 0
    ψ[end] =1
end
    res = sum(abs.(residuals).^2)
    return y,ψ,res
end
