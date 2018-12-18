#make true gpe
using Pkg, Test
pkg"activate ."
using VortexDistributions, Plots, LaTeXStrings, Revise

Κ = 1
Λ = 0.8248
r = linspace(0,6,100)
y,ψ,res = gpecore(Κ)
plt = plot()
xlims!(0,10*Κ);ylims!(0,1.2)
plot!(y,ψ,label=L"k=1",marker="o")
xlabel!(L"$r/\xi$")
ylabel!(L"$\chi(r/\xi)$")
plot!(y,y*Λ,label=L"\Lambda r/\xi")

for j = 2:10
    y,ψ,res = gpecore(j)
    plot!(y,ψ,label=latexstring("k=$(j)"),marker="o",legend = :bottomright)
    display(plt)
end
plt


y,ψ,res = gpecore(Κ)
lt = plot()
plot!(y,ψ,label=L"k=1",marker="o")
xlims!(0,10*Κ)

ψi = make_fastcore(2)

x = linspace(0,10,40)
ψint = ψi(x)

plot!(x,ψint)

ψ2 = vortexcore(x,ψi)
plot!(x,ψ2)

@test ψ2 == ψint
