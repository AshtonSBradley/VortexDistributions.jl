function circmask(ψ,x,y,R)
    for j in eachindex(x), k in eachindex(y)
            (x[j]^2+y[k]^2 < R^2) || (ψ[j,k]=complex(0.))
    end
    return ψ
end
