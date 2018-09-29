function makevortex(ψ,vortex,x,y,ξ=1.0)
    @assert typeof(x)==Array{Float64,1}
    @assert typeof(y)==Array{Float64,1}
    @assert typeof(ψ)==Array{Complex{Float64},2}
    x0, y0, σ0 = vortex
    R(x,y) = sqrt(x^2+y^2)
    return  @. ψ*vortexcore(R(x.-x0,y'.-y0),ξ)*exp(im*σ0*atan(x.-x0,y'.-y0))
end

function makevortex!(ψ,vortex,x,y,ξ=1.0)
    @assert typeof(x)==Array{Float64,1}
    @assert typeof(y)==Array{Float64,1}
    @assert typeof(ψ)==Array{Complex{Float64},2}
    x0, y0, σ0 = vortex
    R(x,y) = sqrt(x^2+y^2)
    ψ .= @. ψ*vortexcore(R(x.-x0,y'.-y0),ξ)*exp(im*σ0*atan(x.-y0,y'.-x0))
end
