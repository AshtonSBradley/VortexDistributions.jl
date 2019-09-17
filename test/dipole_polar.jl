# doubly periodic dipole
# https://dx.doi.org/10.1103/PhysRevLett.112.145301
# should be eigenstate of Hamiltonian evolution

using Test, Plots, Revise, VortexDistributions
gr()

# Methods
H(x) = x > 0. ? 1.0 : 0.0
T(x,xi) = x - xi
tans(x,xk) = tan((T(x,xk) - π)*0.5)
tanhs(x,xk,j) = tanh((T(x,xk) + 2*π*j)*0.5)

function Krow(x,y,xp,yp,xn,yn,j)
    return atan(tanhs(y,yn,j)*tans(x,xn)) -
    atan(tanhs(y,yp,j)*tans(x,xp))
end

# Dimensionless form
function θd(x,y,dip)
    vp,vn = dip
    xp,yp,_ = rawData(vp)
    xn,yn,_ = rawData(vn)
    s = 0.0
    for j = -5:5
        s += Krow(x,y,xp,yp,xn,yn,j)
    end
    return s + π*(H(T(x,xp)) - H(T(x,xn))) - y*(xp - xn)/(2*π)
end

# Now for arbitrary domains and dipole sizes:
function thetad(x,y,xp,yp,xn,yn)
    s = 0.0
    for j = -5:5
        s += Krow(x,y,xp,yp,xn,yn,j)
    end
    s += π*(H(T(x,xp)) - H(T(x,xn))) - y*(xp - xn)/(2*π)
    return s - x*H(abs(yp - yn) - π) + y*H(abs(xp - xn) - π)
end

function Thetad(x,y,xp,yp,xn,yn)
    Lx = x[end]-x[1]
    Ly = y[end]-y[1]
    return @. angle(exp(im*thetad.(x*2*pi/Lx,y'*2*pi/Ly,xp*2*pi/Lx,yp*2*pi/Ly,xn*2*pi/Lx,yn*2*pi/Ly)))
end

# Sukla's parameters
g=0.1
μ=15.
n0 = μ/g
γ=0.01
x = LinRange(-50,50,1024)
y = LinRange(-50,50,1024)

# dipole
xp = 0.
yp = 5

xn = 0.
yn = -5

vp = PointVortex(xp,yp,1)
vn = PointVortex(xn,yn,-1)
ξ = sqrt(1/μ)
dip = ScalarVortex(ξ,[vp;vn])

phase_periodic = Thetad(x,y,xp,yp,xn,yn)
heatmap(x,y,phase_periodic,transpose=true)

# make dipole wavefunction
psi = Torus(sqrt(n0)*one.(x*y'),x,y)
vortex!(psi,dip)

psip = @. abs(psi.ψ)*exp(im*phase_periodic)

heatmap(x,y,abs2.(psip),transpose=true)
heatmap(x,y,angle.(psip),transpose=true)

using FourierGPE

N = (1024,1024)
L = (100,100)
sim = Sim(L,N)
@unpack X,K = sim

phi = kspace(psip,sim) |> fftshift

DX,DK = dfftall(X,K)
kx,ky = K .|> fftshift

data = @. log(abs2(phi)+eps())
heatmap(kx,ky,data)
plot(data[512,:])

# define polar coordinates
Nkx = length(kx)
Nk = Nkx/2 |> Int
Nθk = 2*Nkx
k = LinRange(0,last(kx),Nk)'
θk = LinRange(0,2*pi,Nθk)

function cart2pol(psi,x)
    Nx = length(x)
    Nr = Nx/2 |> Int
    Nθ = 2*Nx
    r = LinRange(0,last(x),Nr)'
    θ = LinRange(0,2*pi,Nθ)

    psiF = fft(psi)

    # precompute modes
    hx = exp.(im*2*pi*()()/Nx)
    hy = exp.(im*)

    # transform
    psipol = zeros(eltype(psi),Nθ,Nr)
    for m = 1:Nx
        for n = 1:Nx
            Hx = @view hx[:,:,m]
            Hy = @view hy[:,:,n]
            @. psipol += Hx * psiF[m,n] * Hy
        end
    end
    return psipol,θ,r
end

function polar(psi,x,ω,hx,hy)
    N = size(hx)[3]
    p = Progress(N,1,"Polar transform...",20)
    H = Oscillator(N,ω)(x)
    T = inv(H'*H)
    F̂ = T*H'*psi*H*T

    Nx = length(x)
    θ = LinRange(0,2*pi,2*Nx)
    r = LinRange(0,last(x),Nx/2 |> Int)'

    psiFP = zeros(eltype(psi),2*Nx,Nx/2 |> Int)

    for m = 1:N
        next!(p)
        for n = 1:(N + 1 - m)
        Hx = @view hx[:,:,m]
        Hy = @view hy[:,:,n]
        @. psiFP += Hx * F̂[m,n] * Hy
    end
    end
    return psiFP
end

phiP = polar(phiF,kx,ω,hx,hy)

densFP = abs2.(phiP)
heatmap(θk,k',densFP,transpose=true)
xlabel!(L"\theta_k");ylabel!(L"k")

dk = k[2] - k[1]
dθk = θk[2] - θk[1]

densk = sum(densFP;dims=1)'*dθk

plot(k',densk)
