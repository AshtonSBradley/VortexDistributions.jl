# VortexDistributions

[![Build Status](https://travis-ci.org/AshtonSBradley/VortexDistributions.jl.svg?branch=master)](https://travis-ci.org/AshtonSBradley/VortexDistributions.jl)  [![Coverage Status](https://coveralls.io/repos/AshtonSBradley/VortexDistributions.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/AshtonSBradley/VortexDistributions.jl?branch=master)  [![codecov.io](http://codecov.io/github/AshtonSBradley/VortexDistributions.jl/coverage.svg?branch=master)](http://codecov.io/github/AshtonSBradley/VortexDistributions.jl?branch=master)

Tools for working with distributions of two-dimensional quantum vortices in Bose-Einstein condendates.

- [x] Fast, accurate vortex detection.
  - Uses a highly optimized version of the plaquette method (phase intergral around each 4-point plaquette), with recursive interpolation to acheive a good balance between speed and accuracy.
  - At present only tests for charge +/-1 in 2D

- [x] Vortex creation
  - Solves the 2D GPE problem for charge n on the infinite domain
  - Interpolates vortex solution to density and phase imprint on arbitrary 2D domains
- [ ] Compressible/incompressible decomposition
- [ ] Recursive cluster algorithm
- [ ] Vortex correlation functions

# Detection Example
```julia
using VortexDistributions, Plots
gr(transpose=true,xlabel="x",ylabel="y",legend=false)

Lx=200;Nx=400;
Ly=200;Ny=400
x = linspace(-Lx/2,Lx/2,Nx);dx=diff(x)[1]
y = linspace(-Ly/2,Ly/2,Ny);dy=diff(y)[1]

#make vortex near the orgin
facx,facy = rand(2)
testvort = [10+dx*facx 10+dy*facy 1.0]

#construct vortex wavefunction
ψ = one.(x.*y') |> complex
makevortex!(ψ,testvort,x,y);
```

Find all vortices (removing edge vortices by default)
```julia
nt,np,nn,vortices = findvortices(ψ,x,y)
```

Plot phase at successive zoom levels with vortex location and detected location:

```julia
p1=heatmap(x,y,angle.(ψ))
scatter!([vortices[1]], [vortices[2]])
scatter!([testvort[1]], [testvort[2]],m=:cross,ms=10,c=:white)
zr = 215:225
p2=heatmap(x[zr],y[zr],angle.(ψ[zr,zr]))
scatter!([vortices[1]], [vortices[2]])
scatter!([testvort[1]], [testvort[2]],m=:cross,ms=10,c=:white)
zr = 221:221
p3=heatmap(x[zr],y[zr],angle.(ψ[zr,zr]),transpose=true)
scatter!([vortices[1]], [vortices[2]],legend=false)
scatter!([testvort[1]], [testvort[2]],m=:cross,ms=10,c=:white)
p=plot(p1,p2,p3,layout=(1,3),size=(700,200))
```
![](/examples/phase.png)

Plot density at successive zoom levels with vortex location and detected location:

```julia
q1=heatmap(x,y,abs2.(ψ),c=:viridis)
scatter!([vortices[1]], [vortices[2]])
scatter!([testvort[1]], [testvort[2]],m=:cross,ms=10,c=:white)
zr = 210:230
q2=heatmap(x[zr],y[zr],abs2.(ψ[zr,zr]),c=:viridis)
scatter!([vortices[1]], [vortices[2]])
scatter!([testvort[1]], [testvort[2]],m=:cross,ms=10,c=:white)
zr = 221:221
q3=heatmap(x[zr],y[zr],angle.(ψ[zr,zr]),c=:viridis)
scatter!([vortices[1]], [vortices[2]],legend=false)
scatter!([testvort[1]], [testvort[2]],m=:cross,ms=10,c=:white)
q=plot(q1,q2,q3,layout=(1,3),size=(700,200))
```
![](/examples/phase.png)

In this example the vortex was created at
```julia
julia> testvort
Out[25]:
1×3 Array{Float64,2}:
 10.3234  10.314  1.0
 ```
 and found at
 ```julia
 julia> vortices
Out[26]:
1×3 Array{Float64,2}:
 10.3274  10.3178  1.0
 ```

 The benchmark gives (2015 MacBook Air 1.6GHz i5)
 ```julia
 using BenchmarkTools
 julia> @btime nt,np,nn,vortices = findvortices(ψ,x,y)
  6.165 ms (550 allocations: 3.96 MiB)
 ```

#### Acknowledgements
Matthew Reeves and Thomas Billam are the primary co-authors of this package.

# External links
___Signatures of Coherent Vortex Structures in a Disordered 2D Quantum Fluid___,\
Matthew T. Reeves, Thomas P. Billam, Brian P. Anderson, and Ashton S. Bradley, \
[Physical Review A 89, 053631 (2014)](http://journals.aps.org/pra/abstract/10.1103/PhysRevA.89.053631)

___Onsager-Kraichnan Condensation in Decaying Two-Dimensional Quantum Turbulence___,\
Thomas P. Billam, Matthew T. Reeves, Brian P. Anderson, and Ashton S. Bradley, \
[Physical Review Letters 112, 145301 (2014)](http://dx.doi.org/10.1103/PhysRevLett.112.145301)
