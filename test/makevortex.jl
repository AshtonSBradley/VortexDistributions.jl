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
