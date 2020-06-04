"""
    f(x)=scalar_ansatz(x)

Evaluate the simple vortex core ansatz at radial point `x`.

The ansatz is the rational function approximation to the vortex core
```math
f(x)=\\sqrt{\\frac{x^2}{1+x^2}}
```
"""
scalar_ansatz(x) = sqrt(x^2/(1+x^2))

"""
    core = Ansatz(ψ,ξ,Λ)

Construct a fast interpolation for the vortex core ansatz.

Returns a callable type that can be evaluated at position `x`.

# Arguments
- `ψ::Interpolation=ψa`: ansatz wavefunction.
- `ξ::Float64=1.0`: healing length.
- `Λ::Float64=0.8249`: slope of wavefunction at vortex core.

# Examples
```jldoctest
julia> f = Ansatz()
julia> x = LinRange(0,10,100)
julia> y = f.(x)
```

See also: [`Exact`](@ref)
"""
Ansatz() = Ansatz(ψa,1.0,Λ) #offending line precomp

function (core::Ansatz)(x)
    @unpack f,ξ,Λ = core
    return f(Λ*x/ξ)
end
(core::Ansatz)(x,y) = core(hypot(x,y))

"""
    core = Exact(ψ,ξ)

Construct fast interpolation for exact vortex core.

Returns a callable type that can be evaluated at position `x`.

# Examples
```jldoctest
julia> f = Exact()
julia> x = LinRange(0,10,100)
julia> y = f.(x)
```

See also: [`Ansatz`](@ref)
"""
Exact(ξ::Float64) = Exact(VortexDistributions.ψi,ξ::Float64)
Exact() = Exact(1.0) # offending line precomp
function (core::Exact)(x)
    @unpack ξ,f = core
    return f(x/ξ)
end
(core::Exact)(x,y) = core(hypot(x,y))

"""
    vort = ScalarVortex(ξ::Float64=1.0;pv::PointVortex)
    vort = ScalarVortex(pv::Array{PointVortex,1})
    vort = ScalarVortex(pv::PointVortex)
    vort = ScalarVortex(ξ::Float64,pv::Array{PointVortex,1})
    vort = ScalarVortex(ξ::Array{Float64,1},pv::Array{PointVortex,1})

Construct a scalar vortex with healing length `ξ` and point vortex coordinates `pv`.

# Arguments
- `ξ::Float64=1.0`: healing length.
- `pv::PointVortex`: scalar or vector of point vortices.

# Examples
```jldoctest
julia> pv = PointVortex(-2.,3.,1)
julia> v1 = ScalarVortex(pv)
```
returns a scalarvortex `v1` suitable for vortex construction in a wavefunction.

See also: [`vortex!`](@ref), [`rand_scalarvortex`](@ref), [`rand_vortex`](@ref), [`PointVortex`](@ref), [`rand_pointvortex`](@ref)
"""
ScalarVortex(vort::PointVortex) = ScalarVortex(Exact(),vort)
function (s::ScalarVortex{T})(x,y) where T <: CoreShape
    xv,yv,qv = s.vort
    return s.core(x - xv, y - yv)*exp(im*qv*atan(y - yv,x - xv))
end
ScalarVortex(ξ::Float64,pv::Array{PointVortex,1}) = ScalarVortex.([Exact(ξ)],pv::Array{PointVortex,1})
ScalarVortex(ξ::Array{Float64,1},pv::Array{PointVortex,1}) = @. ScalarVortex(Exact(ξ),pv::Array{PointVortex,1})
ScalarVortex(ξ::Float64,pv::PointVortex) = ScalarVortex(ξ,[pv])
ScalarVortex(pv::Array{PointVortex,1}) = ScalarVortex(1.0,pv)

"""
    rv = rand_scalarvortex()
    rv = rand_scalarvortex(n::Int)
    rv = rand_scalarvortex(ψ::Field)
    rv = rand_scalarvortex(n::Int,ψ::Field)

Sample `n` random scalar vortices, using the sptial domain of the field `ψ`.

See also: [`Field`](@ref), [`rand_vortex`](@ref), [`ScalarVortex`](@ref), [`randPointVortex`](@ref)
"""
rand_scalarvortex() = ScalarVortex(rand_pointvortex())
rand_scalarvortex(n) = ScalarVortex.(rand_pointvortex(n))
rand_scalarvortex(psi::Field) = ScalarVortex(Exact(),rand_pointvortex(psi))
rand_scalarvortex(n,psi::Field) = ScalarVortex.([Exact()],rand_pointvortex(n,psi))

rand_vortex() = rand_scalarvortex()
rand_vortex(n) = rand_scalarvortex(n)
rand_vortex(n,psi::Field) = rand_scalarvortex(n,psi)

"""
    vortex!(ψ<:Field, vort<:Vortex)

Density and phase imprint a vortex onto the field `ψ`, writing in place.

See also: [`Field`](@ref), [`ScalarVortex`](@ref), [`PointVortex`](@ref), [`rand_pointvortex`](@ref), [`rand_scalarvortex`](@ref)
"""
function vortex!(psi::F,vort::ScalarVortex{T}) where {T <: CoreShape, F<:Field}
    @unpack ψ,x,y = psi
    @. ψ *= vort(x,y')
    @pack! psi = ψ
end
function vortex!(psi::F,vort::S) where {F <: Field, S <: Vortex}
    vortex!(psi,ScalarVortex(vort))
end
function vortex!(psi::F,vort::Array{S}) where {F <: Field, S <: Vortex}
    for j in eachindex(vort)
        vortex!(psi,vort[j])
    end
end

function rand_vortexfield(n)
    Nx = 400; Ny = 400
    Lx = 200; Ly = 200
    x = LinRange(-Lx / 2, Ly / 2, Nx)
    y = LinRange(-Ly / 2, Ly / 2, Ny)

    psi0 = one.(x*y') |> complex; psi = Torus(psi0,x,y)
    vort = rand_vortex(n,psi)
    vortex!(psi,vort)
    return psi,PointVortex(vort)
end

#--- periodic dipole phase
# Billam et al, PRL 112, 145301 (2014), Supplemental
H(x) = x > 0. ? 1.0 : 0.0
shift(x,xi) = x - xi
tans(x,xk) = tan((shift(x,xk) - π)*0.5)
tanhs(x,xk,j) = tanh((shift(x,xk) + 2*π*j)*0.5)

function kernel(x,y,xp,yp,xn,yn,j)
    return atan(tanhs(y,yn,j)*tans(x,xn)) -
    atan(tanhs(y,yp,j)*tans(x,xp))
end

# arbitrary domains and dipole sizes:
function thetad(x,y,xp,yp,xn,yn)
    s = 0.0
    for j = -5:5
        s += kernel(x,y,xp,yp,xn,yn,j)
    end
    s += π*(H(shift(x,xp)) - H(shift(x,xn))) - y*(xp - xn)/(2*π)
    return s - x*H(abs(yp - yn) - π) + y*H(abs(xp - xn) - π)
end

function dipole_phase(x,y,xp,yp,xn,yn)
    Lx = x[end]-x[1]
    Ly = y[end]-y[1]
    return @. thetad(x*2*pi/Lx,y'*2*pi/Ly,xp*2*pi/Lx,yp*2*pi/Ly,xn*2*pi/Lx,yn*2*pi/Ly)
end

function periodic_dipole!(psi::F,dip::Array{ScalarVortex{T},1}) where {T <: CoreShape, F<:Field}
    @assert length(dip) == 2
    @assert dip[1].vort.q + dip[2].vort.q == 0
    @unpack ψ,x,y = psi
    (dip[1].vort.q > 0) ? (jp = 1;jn = 2) : (jp = 2;jn = 1)
    vp = vortex_array(dip[jp].vort)[1:2]
    vn = vortex_array(dip[jn].vort)[1:2]
    ψ .*= abs.(dip[jn].core.(x,y')*dip[jp].core.(x,y'))
    ψ .*= exp.(im*dipole_phase(x,y,vp...,vn...))
    @pack! psi = ψ
end

function periodic_dipole!(psi::F,dip::Dipole) where F <: Field
    @unpack ψ,x,y = psi
    rp = vortex_array(dip.vp)[1:2]
    rn = vortex_array(dip.vn)[1:2]
    ψ .*= abs.(dip[jn].core.(x,y')*dip[jp].core.(x,y'))
    ψ .*= exp.(im*dipole_phase(x,y,rp...,rn...))
    @pack! psi = ψ
end


#Chebyshev methods

function gpecore_exact(K,L=2,N=100,R = K)
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
y =  @. R*(1+z)/(1-z)

#initial guess based on ansatz
ψ = @. y/hypot(1/Λ,y)
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

    Jacobi = -0.5*( (Qmat.*(D2z) + Zmat.*(Dz) ).*(Qmat/(4*R^2))
    + (Qmat/(2*R).*Ymat).*(Dz) ) + diagm(0 => 0.5*K^2 ./y.^2)+ diagm(0 => 3*ψ.^2 .- 1)


    Jacobi[1,:] = [1 zeros(1,N-1)]
    Jacobi[N,:] = [zeros(1,N-1) 1]


    Δ = vec(-4/7*(Jacobi\residuals))

    ψ = ψ + Δ
    ψ[1] = 0
    ψ[end] =1
end
    res = norm(residuals)^2
    return y,ψ,res
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
