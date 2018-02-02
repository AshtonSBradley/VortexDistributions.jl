function makevortex(ψ,vortex,x,y,ξ)
    x0, y0, σ0 = vortex
    X = x-x0; Y = y'-y0
    R = @. sqrt(X^2+Y^2)
    return  ψ.*vortexcore(R,ξ).*exp.(im*σ0*atan2.(X,Y))
end
