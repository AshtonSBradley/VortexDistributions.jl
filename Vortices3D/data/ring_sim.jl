using FourierGPE

##Set simulation parameters
    L=(16.,16.,16.);
    N=(64,64,64);
    sim = Sim(L,N);
    @unpack_Sim sim;

## Initialse sim
    # parameters
    μ = 25.0;
    γ = 0.05;
    tf = 4/γ;
    Nt = 200;
    t = LinRange(0.,tf,Nt);

## Run sim
    x,y,z = X;
    ψi = randn(N)+im*randn(N);
    ϕi = kspace(ψi,sim);

    @pack_Sim! sim;

## Evolve in k-space
    sol = runsim(sim); # will take a few minutes to run.

# size(xspace(sol[40], sim))

plot_iso(xspace(sol[150], sim))
psi_chaos1 = xspace(sol[10], sim)
psi_chaos2 = xspace(sol[16], sim)
psi_knots1 = xspace(sol[20], sim)
psi_knots2 = xspace(sol[25], sim)
psi_knots3 = xspace(sol[30], sim)
psi_knots4 = xspace(sol[40], sim)
psi_tubes1 = xspace(sol[50], sim)
psi_tubes2 = xspace(sol[90], sim)
psi_ringtube = xspace(sol[100], sim)

@save "box_vorts.jld2" psi_chaos1 psi_chaos2 psi_knots1 psi_knots2 psi_knots3 psi_knots4 psi_tubes1 psi_tubes2 psi_ringtube X

# function dense(i)
#     ψm = xspace(sol[i],sim)
#     density = abs2.(ψm)
#     pmax = maximum(density)
#     return density/pmax
# end

# function densityfilm(Nt,saveto="media/3dquenchiso.gif")
#     scene = Scene()
#     tindex = Node(1)

#     scene = volume(lift(i -> dense(i), tindex), algorithm = :iso,show_axis=false)

#     record(scene, saveto, 1:Nt-10) do i
#         tindex[] = i
#     end
# end

# densityfilm(Nt)