function [a, e, i, OM, om, th] = car2kep(rr, vv, mu)
% 
% Function to convert Cartesian coordinates into Keplerian elements.
% 
% PROTOTYPE:
%  [a, e, i, OM, om, th] = car2kep(rr, vv, mu)
% 
% INPUT:
%  rr [3,1]       Cartesian position vector [km]
%  vv [3,1]       Cartesian velocity vector [km/s]
%  mu [1,1]       Gravitational parameter of the central body [km^3/s^2]
% 
% OUTPUT:
%  a [1,1]        Semi-major axis [km]
%  e [1,1]        Eccentricity [dimensionless]
%  i [1,1]        Inclination [rad]
%  OM [1,1]       Right Ascension of Ascending Node (RAAN) [rad]
%  om [1,1]       Argument of pericenter [rad]
%  th [1,1]       True anomaly [rad]
% 
% AUTHORS: Richero Giovanni, Emrem Mert, Bernarsconi Ludovico, Serlini
% Mariagiulia
%
% -------------------------------------------------------------------------

% Compute the magnitude of the position vector
r = norm(rr);

% Compute the magnitude of the velocity vector
v = norm(vv);

% Compute the semi-major axis using the vis-viva equation
a = (2/r - v^2/mu)^(-1);

% Compute the specific angular momentum vector and its magnitude
hh = cross(rr, vv);
h = norm(hh);

% Compute the eccentricity vector and its magnitude
ee = cross(vv, hh)./mu - rr./r;
e = norm(ee);

% Compute the inclination
i = acos(hh(3)/h);

% Compute the line of nodes
N = cross([0 0 1]', hh) / norm(cross([0 0 1]', hh));

% Handle the special case of zero inclination
if i == 0
    N = [1 0 0]';
end

% Compute the Right Ascension of Ascending Node (RAAN)
OM = acos(N(1)) + (N(2) < 0) * (2*pi - 2*acos(N(1)));

% Compute the Argument of Pericenter
om = acos(N' * ee / e) + (ee(3) < 0) * (2*pi - 2*acos(N' * ee / e));

% Compute the True Anomaly
th = acos(rr' * ee / (r * e)) + (vv' * rr < 0) * (2*pi - 2*acos(rr' * ee / (r * e)));

end


