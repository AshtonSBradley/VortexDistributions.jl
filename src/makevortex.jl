function makevortex(ψ,vortex,x,y,ξ)
    x0, y0, σ0 = vortex
    ϕ = angle.(ψ)
    U = abs.(ψ)
    X = x-x0
    Y = y-y0
    R = @. sqrt(X^2+Y^2)
    ϕ .+= @. σ0*atan2(X,Y)
    χ = @. U*vortexcore(R,ξ)
    return  @. U*exp(im*ϕ)
end
