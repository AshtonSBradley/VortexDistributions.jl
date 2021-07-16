using JLD2, Interpolations
@load joinpath(@__DIR__,"src/exactcore.jld2") ψi
@load joinpath(@__DIR__,"src/ansatzcore.jld2") ψa

psii_coefs = ψi.coefs
psii_knots = ψi.knots
psia_coefs = ψa.coefs
psia_knots = ψa.knots
@save joinpath(@__DIR__,"src/data.jld2") psii_coefs psii_knots psia_coefs psia_knots

## save in new format
@load joinpath(@__DIR__,"src/data.jld2") psii_coefs psii_knots psia_coefs psia_knots
ψi = interpolate(psii_knots, psii_coefs, Gridded(Linear()))
ψa = interpolate(psia_knots, psia_coefs, Gridded(Linear()))
@save joinpath(@__DIR__,"src/cores.jld2") ψi ψa

## load interp and check type 
@load joinpath(@__DIR__,"src/cores.jld2") ψi ψa
