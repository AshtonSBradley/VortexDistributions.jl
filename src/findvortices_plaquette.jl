"""
    vortices = findvortices(x,y,ψ)

 Locates vortices as 2π phase windings around plaquettes on a cartesian spatial field.

Requires a 2D wavefunction ψ. Non-optimized version encountered in previous literature.
"""

function findvortices_plaquette(ψ,x,y)
@assert typeof(x)==Array{Float64,1}
@assert typeof(y)==Array{Float64,1}
@assert typeof(ψ)==Array{Complex{Float64},2}

phase = angle.(ψ)
Nx,Ny = size(ψ)

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

        ixp,iyp,vp = findnz(vortexgrid.>0.)
        xp = x[ixp]; yp = y[iyp]

        ixn,iyn,vn = findnz(vortexgrid.<0.)
        xn = x[ixn]; yn = y[iyn];

        vortices = [xn yn -vn; xp yp vp] |> sortrows
        return vortices
end
