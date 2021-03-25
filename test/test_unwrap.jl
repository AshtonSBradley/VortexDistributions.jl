##
using Plots
# import VortexDistributions:unwrap
using VortexDistributions

##
# function unwrap(phase::Array{Float64,1})
#     uphase = copy(phase)
#     s1 = length(phase)
#         @inbounds (k=uphase[1] - uphase[s1];k > π) && (uphase[1] -= 2π*div(k,2pi,RoundNearest))
#         @inbounds (k=uphase[1] - uphase[s1];k < -π) && (uphase[1] += 2π*div(-k,2pi,RoundNearest))
#         for i in 2:s1
#         @inbounds (k=uphase[i] - uphase[i-1];k > π) && (uphase[i] -= 2π*div(k,2pi,RoundNearest))
#         @inbounds (k=uphase[i] - uphase[i-1];k < -π) && (uphase[i] += 2π*div(-k,2pi,RoundNearest))
#         end
#         return uphase
# end

##
N=1000
θ = LinRange(0,4*pi,N)
psi = exp.(im*θ) + 0.1*(randn(N)+im*randn(N))
ϕ = angle.(psi)
ϕu = unwrap(ϕ)
ϕui = similar(ϕu)
unwrap!(ϕui,ϕ)

#
plot(ϕ,label="phi",c=:red,alpha=0.4,lw=6)
plot!(ϕu,label="phiu",lw=2,c=:black)
plot!(ϕui,label="phiu",lw=2,c=:green)