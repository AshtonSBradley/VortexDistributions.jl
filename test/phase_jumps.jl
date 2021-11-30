using VortexDistributions

# simple test field
Nx = 400; Ny = 400
Lx = 200; Ly = Lx
x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, Ny+1)[1:end-1];
psi0 = one.(x*y') |> complex
psi = Torus(copy(psi0),x,y)
xp = 1.133
yp = 2.787
vp = PointVortex(xp,yp,1)
sp = ScalarVortex(vp)
vortex!(psi,sp)

ψ = psi.ψ
phase = angle.(ψ)

diffx,diffy = phase_jumps(phase, 1),phase_jumps(phase, 2) 
a = zero(phase)
b = zero(phase) 
phase_jumps!(a, phase)
@test a == diffx
phase_jumps!(b, phase, 2)
@test b == diffy

