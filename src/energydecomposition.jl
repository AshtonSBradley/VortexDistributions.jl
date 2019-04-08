function energydecomp(ψ,kx,ky,kz,k2)
    rho = abs2.(ψ)
    vx,vy,vz = velocities(ψ,kx,ky,kz,k2)
    wx = @. sqrt(rho)*vx; wy = @. sqrt(rho)*vy; wz = @. sqrt(rho)*vx
    wxi,wyi,wzi,wxc,wyc,wzc = helmholtzdecomp(wx,wy,wz,kx,ky,kz,k2)
    et = @. abs2(wx)+abs2(wy)+abs2(wz); et *= 0.5
    ei = @. abs2(wxi)+abs2(wyi)+abs2(wzi); ei *= 0.5
    ec = @. abs2(wxc)+abs2(wyc)+abs2(wzc); ec *= 0.5
    return et,ei,ec
end

function energydecomp(ψ,kx,ky,k2)
    rho = abs2.(ψ)
    vx, vy = velocities(ψ,kx,ky,k2)
    wx = @. sqrt(rho)*vx; wy = @. sqrt(rho)*vy
    wxi, wyi, wxc, wyc = helmholtzdecomp(wx,wy,kx,ky,k2)
    et = @. abs2(wx)+abs2(wy); et *= 0.5
    ei = @. abs2(wxi)+abs2(wyi); ei *= 0.5
    ec = @. abs2(wxc)+abs2(wyc); ec *= 0.5
    return et,ei,ec
end
