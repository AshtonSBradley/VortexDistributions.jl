# VortexDistributions

[![Build Status](https://travis-ci.org/AshtonSBradley/VortexDistributions.jl.svg?branch=master)](https://travis-ci.org/AshtonSBradley/VortexDistributions.jl)  [![Coverage Status](https://coveralls.io/repos/AshtonSBradley/VortexDistributions.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/AshtonSBradley/VortexDistributions.jl?branch=master)  [![codecov.io](http://codecov.io/github/AshtonSBradley/VortexDistributions.jl/coverage.svg?branch=master)](http://codecov.io/github/AshtonSBradley/VortexDistributions.jl?branch=master)

Tools for working with distributions of two-dimensional quantum vortices in Bose-Einstein condendates.

- [x] Fast, accurate vortex detection.
  - Uses a highly optimized version of the plaquette method (phase intergral around each 4-point plaquette), with recursive interpolation to acheive a good balance between speed and accuracy. 
  - At present only tests for charge +/-1 in 2D
- [ ] Vortex creation
  - Solves the 2D GPE problem for charge n on the infinite domain
  - Interpolates vortex solution to density and phase imprint on arbitrary 2D domains
- [ ] Compressible/incompressible decomposition
- [ ] Recursive cluster algorithm
- [ ] Vortex correlation functions

#### Acknowledgements
Matthew Reeves provided the code for solving the charge n vortex problem.

# External links
___Signatures of Coherent Vortex Structures in a Disordered 2D Quantum Fluid___,\
Matthew T. Reeves, Thomas P. Billam, Brian P. Anderson, and Ashton S. Bradley, \
[Physical Review A 89, 053631 (2014)](http://journals.aps.org/pra/abstract/10.1103/PhysRevA.89.053631)

___Onsager-Kraichnan Condensation in Decaying Two-Dimensional Quantum Turbulence___,\
Thomas P. Billam, Matthew T. Reeves, Brian P. Anderson, and Ashton S. Bradley, \
[Physical Review Letters 112, 145301 (2014)](http://dx.doi.org/10.1103/PhysRevLett.112.145301)
