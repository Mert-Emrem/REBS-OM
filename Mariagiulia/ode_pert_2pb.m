function dy = ode_pert_2pb(~, y, mu, R_E, J_2)
%ode_2bp ODE system for the two-body problem (Keplerian motion)
%
% PROTOTYPE
% dy = ode_2bp( t, y, mu )
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% y[6x1] State of the body ( rx, ry, rz, vx, vy, vz ) [ L, L/T ]
% mu[1] Gravitational parameter of the primary [L^3/T^2]
%
% OUTPUT:
% dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]

r = y(1:3);
v = y(4:6);
rnorm = norm(r);

a = 3/2 * J_2 * R_E^2/rnorm^4 * [r(1)/rnorm*(5*r(3)/rnorm-1); 
                                 r(2)/rnorm*(5*r(3)/rnorm-1);
                                 r(3)/rnorm*(5*r(3)/rnorm-3)];

dy = [v(1); v(2); v(2); (-mu/rnorm^3)*r(1) + a(1); (-mu/rnorm^3)*r(2) + a(2); (-mu/rnorm^3)*r(1) + a(3)];

end