#make true gpe

using VortexDistributions, Plots, LaTeXStrings

Κ=1
Λ = 0.8248
r = linspace(0,6,100)
y,ψ,res = gpecore(r,Κ)
plt = plot()
plot!(y,ψ,label=L"k=1",marker="o")
xlims!(0,10*Κ)
xlabel!(L"$r/\xi$")
ylabel!(L"$\chi(r/\xi)$")
#plot!(y,abs2.(ψ),label=L"|\chi(r/\xi)|^2")
plot!(y,y*Λ,label=L"\Lambda r/\xi")
ylims!(0,1.2)
for j = 2:10
    y,ψ,res = gpecore(r,1.0,j)
    plot!(y,ψ,label=latexstring("k=$(j)"),marker="o",legend = :bottomright)
    display(plt)
end
plt


y,ψ,res = gpecore(r,1.0,Κ,2,30,Κ)
lt = plot()
plot!(y,ψ,label=L"k=1",marker="o")
xlims!(0,10*Κ)
