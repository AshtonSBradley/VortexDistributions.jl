scalaransatz(x) = sqrt(x^2/(1+x^2))
r(x,y) = sqrt(x^2+y^2)

Ansatz() = Ansatz(ψa,1.0,Λ)

function (core::Ansatz)(x)
    @unpack f,ξ,Λ = core
    return f(Λ*x/ξ)
end
(core::Ansatz)(x,y) = core(r(x,y))

Exact() = Exact(ψi,1.0)

function (core::Exact)(x)
    @unpack ξ,f = core
    return f(x/ξ)
end
(core::Exact)(x,y) = core(r(x,y))

ScalarVortex(vort::PointVortex) = ScalarVortex(Exact(),vort)

function (s::ScalarVortex{T})(x,y) where T <: CoreShape
    @unpack xv,yv,qv = s.vort
    return s.core(x - xv, y - yv)*exp(im*qv*atan(y - yv,x - xv))
end

randScalarVortex() = ScalarVortex(randPointVortex())
randScalarVortex(n) = ScalarVortex.(randPointVortex(n))
randScalarVortex(psi::Field) = ScalarVortex(Exact(),randPointVortex(psi))
randScalarVortex(n,psi::Field) = ScalarVortex.([Exact()],randPointVortex(n,psi))

randVortex() = randScalarVortex()
randVortex(n) = randScalarVortex(n)
randVortex(n,psi::Field) = randScalarVortex(n,psi)

function vortex!(psi::F,vort::ScalarVortex{T}) where {T <: CoreShape, F<:Field}
    @unpack ψ,x,y = psi
    @. ψ *= vort.(x,y')
    @pack! psi = ψ
end

function vortex!(psi::F,vort::Array{S}) where {F <: Field, S <: Vortex}
    for j in eachindex(vort)
        vortex!(psi,vort[j])
    end
end

function randVortexField(n)
    Nx = 400; Ny = 400
    Lx = 200; Ly = 200
    x = LinRange(-Lx / 2, Ly / 2, Nx)
    y = LinRange(-Ly / 2, Ly / 2, Ny)

    psi0 = one.(x*y') |> complex; psi = Torus(psi0,x,y)
    vort = randVortex(n,psi)
    vortex!(psi,vort)
    return psi,PointVortex(vort)
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
y =  R*(1 .+z)./(1 .-z)

#initial guess based on ansatz
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
