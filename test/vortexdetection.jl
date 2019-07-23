# detection

N = 200
L = 100.0
x = LinRange(-L/2,L/2,N)
y = x
psi0 = one.(x*y') |> complex
psi = Torus(psi0,x,y)

vtest = randPointVortex(1,psi)
gpeansatz!(psi,vtest)
@unpack ψ = psi
heatmap(x,y,angle.(ψ),transpose=true)
