function velocities(ψ,kx,ky,k2)
    rho = abs2.(ψ)
    ϕ = fft(ψ)
    ψx = ifft(im*kx.*ϕ)
    ψy = ifft(im*ky.*ϕ)
    vx = @. 0.5*im*(conj(ψx)*ψ-conj(ψ)*ψx)/rho
    vy = @. 0.5*im*(conj(ψy)*ψ-conj(ψ)*ψy)/rho
    return vx, vy
end

function velocities(ψ,kx,ky,kz,k2)
    rho = abs2.(ψ)
    ϕ = fft(ψ)
    ψx = ifft(im*kx.*ϕ)
    ψy = ifft(im*ky.*ϕ)
    ψz = ifft(im*kz.*ϕ)
    vx = @. 0.5*im*(conj(ψx)*ψ-conj(ψ)*ψx)/rho
    vy = @. 0.5*im*(conj(ψy)*ψ-conj(ψ)*ψy)/rho
    vz = @. 0.5*im*(conj(ψz)*ψ-conj(ψ)*ψz)/rho
    return vx,vy,vz
end
