#make true gpe
#using Pkg
#pkg"activate ."


using VortexDistributions, Plots, LaTeXStrings

Κ=2
Λ = 0.8248

y,ψ,res = gpecore(Κ)
plt = plot()
plot!(y,ψ,label=L"k=1")
xlims!(0,10*Κ)
xlabel!(L"$r/\xi$")
ylabel!(L"$\chi(r/\xi)$")
#plot!(y,abs2.(ψ),label=L"|\chi(r/\xi)|^2")
plot!(y,y*Λ,label=L"\Lambda r/\xi")
ylims!(0,1.2)
for j = 2:10
    y,ψ,res = gpecore(j)
    plot!(y,ψ,label=latexstring("k=$(j)"),legend = :bottomright)
    display(plt)
end
plt
