# vortex core construction

function vortexcore(r,ξ,ansatz=true,charge=1)
    ansatz ? (return r/sqrt(r^2 + ξ^2)) : (return gpecore(r/ξ,charge))
end

function gpecore(r,Κ,L=2,N=100,R=Κ)
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
y =  R*(1 .+z)./(1 .-z)  ;

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
        + (Q/(2*R) ./y).*(Dz*ψ) )+ 0.5*Κ^2 *ψ ./y.^2 + ψ.^3 .- ψ
residuals[1] = 0; residuals[end] =0

while sum(abs.(residuals).^2) > 1e-12
    # Second Derivative
    residuals = -0.5*( (Q.*(D2z*ψ) + 2*(z .-1).*(Dz*ψ) ).*Q/(4*R^2)
            + (Q/(2*R) ./y).*(Dz*ψ) )+ 0.5*Κ^2 *ψ ./y.^2 + ψ.^3 .- ψ
    residuals[1] = 0; residuals[end] =0

    #println(sum(abs.(residuals).^2))


    Jacobi = -0.5*( (Qmat.*(D2z) + Zmat.*(Dz) ).*(Qmat/(4*R^2))
    + (Qmat/(2*R).*Ymat).*(Dz) ) + diagm(0 => 0.5*Κ^2 ./y.^2)+ diagm(0 => 3*ψ.^2 .- 1)


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
