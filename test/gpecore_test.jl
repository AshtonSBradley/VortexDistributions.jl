#make true gpe

using VortexDistributions, Plots, LaTeXStrings
gr();
Κ=1
Λ = 0.8248
r = linspace(0,6,100)
y,ψ,res = gpecore(r,Κ)
plt = plot()
xlims!(0,10*Κ);ylims!(0,1.2)
plot!(y,ψ,label=L"k=1",marker="o")
xlabel!(L"$r/\xi$")
ylabel!(L"$\chi(r/\xi)$")
plot!(y,y*Λ,label=L"\Lambda r/\xi")

for j = 2:10
    y,ψ,res = gpecore(r,j)
    plot!(y,ψ,label=latexstring("k=$(j)"),marker="o",legend = :bottomright)
    display(plt)
end
plt


y,ψ,res = gpecore(r,Κ)
lt = plot()
plot!(y,ψ,label=L"k=1",marker="o")
xlims!(0,10*Κ)
