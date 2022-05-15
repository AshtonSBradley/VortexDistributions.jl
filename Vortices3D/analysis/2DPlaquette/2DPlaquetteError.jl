include("util.jl")
include("2Danalysis.jl")

n = 16
Lx = 10; Ly = Lx;
Nx = n; Ny = n;
x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, Ny+1)[1:end-1];
# x = LinRange(-Lx/2,Lx/2, Nx); y = LinRange(-Ly/2,Ly/2, Ny);
psi0 = one.(x*y') |> complex
psi = Torus(psi0,x,y)



v = PointVortex(uniform(x[2], x[end-1]), uniform(y[2],y[end-1]), 1)
v2 = PointVortex(uniform(x[2], x[end-1]), uniform(y[2],y[end-1]), 1)

vortex!(psi, v)
vortex!(psi, v2)


##
vf = findvortices(psi)

findMinAvgError(psi, [v, v2], vf)
findMinAvgError2(psi, [v, v2], vf)
errorTest(2, 128)

@time errvec = [x = errorTest(10, 128, 2) for i in 1:10000 if x!=400]
mean(errvec)
findmax(errvec)

##
vf = findvortices(psi)
naivePlaquette(psi)
remove_vortices_edge(findvortices_jumps(psi), psi)

dx = x[2]-x[1]
euclidVorts(v, vf[1])/dx

maxError = sqrt(2)*(dx/((30/4)))

function runErrorSim(Ntrials)
    diffs = []
    count = 0
    while (count < Ntrials)
        n = 16
        Lx = 10; Ly = Lx;
        Nx = n; Ny = n;
        x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, Ny+1)[1:end-1];
        # x = LinRange(-Lx/2,Lx/2, Nx); y = LinRange(-Ly/2,Ly/2, Ny);
        psi0 = one.(x*y') |> complex
        psi = Torus(psi0,x,y)

        v = PointVortex(uniform(x[2], x[end-1]), uniform(y[2],y[end-1]), 1)
        vortex!(psi, v)
        vf = naivePlaquette(psi.ψ, [x, y], 0)
        push!(diffs, euclidVorts(v, vf[1]))
        count+=1
    end
    return diffs
end

errorBound = sqrt(1/2)*dx
avgError = dx*(1/sqrt(2*pi))


N = 1000000

diffs = runErrorSim(N)
maxError = findmax(diffs)[1]



avgResult = mean(diffs)

maxError/errorBound
avgResult/avgError