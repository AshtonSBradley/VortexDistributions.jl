function helmholtzdecomp(wx,wy,kx,ky,k2)
    wxk = fft(wx); wyk = fft(wy)
    kdotw = @. kx*wxk + ky*wyk
    wxkc = @. kdotw*kx/k2; wxkc[1,1] = zero(wxkc[1,1])
    wykc = @. kdotw*ky/k2; wykc[1,1] = zero(wykc[1,1])
    wxki = @. wxk - wxkc
    wyki = @. wyk - wykc
    wxc = ifft(wxkc); wyc = ifft(wykc)
    wxi = ifft(wxki); wyi = ifft(wyki)
    return wxi, wyi, wxc, wyc
end

function helmholtzdecomp(wx,wy,wz,kx,ky,kz,k2)
    #TODO check k sizes and do some careful tests!
    # ky = ky'
    # kz = reshape(kz,1,1,length(ky))
    wxk = fft(wx); wyk = fft(wy); wzk = fft(wz);
    kdotw = @. kx*wxk + ky*wyk + kz*wzk
    wxkc = @. kdotw*kx/k2; wxkc[1,1] = zero(wxkc[1,1])
    wykc = @. kdotw*ky/k2; wykc[1,1] = zero(wykc[1,1])
    wzkc = @. kdotw*kz/k2; wzkc[1,1] = zero(wzkc[1,1])
    wxki = @. wxk - wxkc
    wyki = @. wyk - wykc
    wzki = @. wzk - wzkc
    wxc = ifft(wxkc); wyc = ifft(wykc); wzc = ifft(wzkc)
    wxi = ifft(wxki); wyi = ifft(wyki); wzi = ifft(wzki)
    return wxi, wyi, wzi, wxc, wyc, wzc
end
