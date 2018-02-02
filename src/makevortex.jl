function makevortex(ψ,vortex,x,y,ξ)
    @assert typeof(x)==Array{Float64,1}
    @assert typeof(y)==Array{Float64,1}
    @assert typeof(ψ)==Array{Complex{Float64},2}
    x0, y0, σ0 = vortex
    X = x-x0; Y = y'-y0
    R = @. sqrt(X^2+Y^2)
    return  ψ.*vortexcore(R,ξ).*exp.(im*σ0*atan2.(Y,X))
end
