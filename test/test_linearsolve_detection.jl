#--- test linear solve method (Foster PhD)
using Plots, Revise, VortexDistributions

pv1 = PointVortex(0.1,0.3,1)
v1 = ScalarVortex(pv1)

x = LinRange(-10,10,256)
y = x

psi = Torus(one.(x.*y') |> complex,x,y)

vortex!(psi,v1)

psiv = psi.ψ

heatmap(x,y,angle.(psiv))

#--- define linear solve (Foster)

function vortex_detect2d(psi)
# Detect vortices in a field psi, on plaquettes of a 2D grid.
#
# First we localize vortices on the grid plaquettes using the phase winding
# algorithm.  After we've found which grid plaquettes contain vortices, we
# estimate the positions more accurately by solving the equations
#
#   real(psi(x,y)) = 0, imag(psi(x,y)) = 0
#
# In principle we could do this using the exact psi as computed via spectral
# interpolation (provided psi is defined in a spectral basis) but that would
# probably be relatively slow.  Instead, we make a linear approximation over
# the plaquette for the real and imaginary parts of psi, and find the position
# where these approximations are mutually zero.  This corresponds to the
# position of the vortex.
#
# The linear approximations are computed from the values of psi at the corners
# of the plaquette using least squares.
#
# Inputs:
#   psi - complex wavefunction in which to detect vortices.
#
# Outputs:
#   positions - Array of vortex positions with shape 2xN where N is the number
#               of vortices:  [ x1  x2  x3 ... xN ]
#                             [ y1  y2  y3 ... yN ]
#   winding - The vortex winding number of the positions array, +1 or -1.

vorticity = phase_winding2d(angle(psi));
[iy,ix] = find(vorticity);
positions = [ix,iy]';

if nargout > 1
    winding = vorticity(vorticity~=0);
end

# % Compute least squares estimator for the linear function coefficients on a
# % plaquette.  Suppose the real part of psi on a plaquette containing a vortex
# % is
# %
# %   [r1 r2;
# %    r3 r4]
# %
# % We want to fit a plane (a*x + b*y + c) to this such that
# %
# %   r1 = a*0 + b*0 + c*1
# %   r2 = a*1 + b*0 + c*1
# %   r3 = a*0 + b*1 + c*1
# %   r4 = a*1 + b*1 + c*1
# %
# % In matrix form, that would be L * [a;b;c] = [r1;r2;r3;r4], which is the
# % origin of the mysterious looking matrix below.
# %
L = [0 0 1;
     1 0 1;
     0 1 1;
     1 1 1];
lsEst = inv(L'*L)*L';
# % (For future reference, lsEst actually turns out to be the following)
# % lsEst = 1/4 * [-2,  2, -2,  2;
# %                -2, -2,  2,  2;
# %                 3,  1,  1, -1];
for i = 1:length(ix)
    # Extract relevant grid plaquette containing the vortex.
    p = psi(iy(i)+[0,1], ix(i)+[0,1]).'
    # Compute linear approximations of real & imag parts of p
    rcoeff = lsEst*real(p[:])
    icoeff = lsEst*imag(p[:])
    # Find position where these are mutually zero
    pos = [rcoeff(1:2), icoeff(1:2)]' \ -[rcoeff(3); icoeff(3)]
    positions[:,i] = positions[:,i] + pos;
    return positions, winding
end
