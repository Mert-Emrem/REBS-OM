function dy = ode_2bp_J2(~, y, MU, J2, Re)
% ODE_2BP_J2 ODE system for the two-body problem with J2 perturbation
%
% This function computes the derivatives of the state vector for the
% two-body problem, including the perturbative effects of the Earth's
% oblateness (J2 effect).
%
% STRUCTURE:
%  dy = ode_2bp_J2(t, y, mu, J2, Re)
%
% INPUT:
%  - t      [1x1]   Time (can be omitted, as the system is autonomous) [T]
%  - y      [6x1]   State vector (rx, ry, rz, vx, vy, vz) [L, L/T]
%  - mu     [1x1]   Gravitational parameter of the primary body [L^3/T^2]
%  - J2     [1x1]   Second zonal harmonic coefficient of the primary body [-]
%  - Re     [1x1]   Equatorial radius of the primary body [L]
%
% OUTPUT:
%  - dy [6x1] Derivative of the state vector [L/T, L/T^2]
%
% CONTRIBUTORS:
%   Ludovico Bernasconi
%
% -------------------------------------------------------------------------

% Position, velocity 
r       = y(1:3);
v       = y(4:6);

% Calculate the second zonal harmonical perturabation a_J2
rnorm   = norm(r);
z       = r(3);
a_J2    = (3/2) * J2 * (MU * Re^2 / rnorm^5) * [(5 * (z/rnorm)^2 - 1) * r(1); (5 * (z/rnorm)^2 - 1) * r(2); (5 * (z/rnorm)^2 - 3) * r(3)];

% Set the derivatives of the state
dy      = [v; (-MU/rnorm^3) * r + a_J2];

end

