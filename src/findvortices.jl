"""
    vortices = findvortices(x,y,ψ)
    
 Locates vortices as 2π phase windings around plaquettes on a cartesian spatial field.

Requires a 2D wavefunction ψ; for 3D, pass slices.
"""

function findvortices(x,y,ψ)
@assert typeof(x)==Array{Float64,1}
@assert typeof(y)==Array{Float64,1}
@assert typeof(ψ)==Array{Complex{Float64},2}

phase = angle.(ψ)
Nx,Ny = size(ψ)
# x corresponds to column of ψ
# y is a row vector
X = x*ones(y')
Y = ones(x)*y'
vortexgrid = zeros(Nx,Ny)

    for i = 1:Nx-1, j = 1:Ny-1
            m = 0;
            Δ = phase[i+1,j]-phase[i,j]
            abs(Δ) > π && (m += sign(-Δ))

            Δ = phase[i+1,j+1]-phase[i+1,j]
            abs(Δ) > π && (m += sign(-Δ))

            Δ = phase[i,j+1]-phase[i+1,j+1]
            abs(Δ) > π && (m += sign(-Δ))

            Δ = phase[i,j]-phase[i,j+1]
            abs(Δ) > π && (m += sign(-Δ))

            vortexgrid[i,j] = m
    end

    pos = vortexgrid .> 0
    neg = vortexgrid .< 0
    indp = find(pos)
    xp = X[indp]; xp = xp[:]
    yp = Y[indp]; yp = yp[:]
    indn = find(neg)
    xn = X[indn]; xn = xn[:]
    yn = Y[indn]; yn = yn[:]
    xv = [xp;xn]; yv=[yp;yn]; s = [ones(xp);-ones(xn)]
    vortices = [xv yv s]
    vortices = sortrows(vortices)
    return vortices
end
