function [rr, vv] = kep2car_vec(a, e, i, OM, om, th, mu)
% Function to convert vectorized Keplerian orbital elements to Cartesian state vectors
% This function calculates the position and velocity vectors in the 
% Earth-Centered Inertial (ECI) frame for multiple sets of Keplerian elements.
% 
% PROTOTYPE:
%  [rr, vv] = kep2car_vec(a, e, i, OM, om, th, mu)
% 
% INPUT:
%  a    [Nx1]  Semi-major axis of the orbit [km]
%  e    [Nx1]  Eccentricity of the orbit [-]
%  i    [Nx1]  Inclination of the orbit [rad]
%  OM   [Nx1]  Right Ascension of the Ascending Node (RAAN) [rad]
%  om   [Nx1]  Argument of perigee [rad]
%  th   [Nx1]  True anomaly [rad]
%  mu   [1x1]  Gravitational parameter of the central body [km^3/s^2]
% 
% OUTPUT:
%  rr   [3xN]  Position vectors in the ECI frame [km]
%  vv   [3xN]  Velocity vectors in the ECI frame [km/s]
% 
% AUTHORS: Richero Giovanni, Emrem Mert, Bernarsconi Ludovico, Serlini Mariagiulia
% -------------------------------------------------------------------------

% Compute semi-latus rectum and specific angular momentum
p = a .* (1 - e.^2);
h = sqrt(p .* mu);

% Calculate the orbital radii
r = p ./ (1 + e .* cos(th));

% Position and velocity in the Perifocal frame
r_PF = r .* [cos(th); sin(th); 0 * sin(th)];
v_PF = (mu ./ h) .* [-sin(th); (e + cos(th)); 0 * sin(th)];

% Initialize rotation matrices for ECI to Perifocal transformation
R_om = zeros(3, 3, length(om));
c_om = cos(om);
s_om = sin(om);
R_om(1, 1, :) = c_om;
R_om(1, 2, :) = s_om;
R_om(2, 1, :) = -s_om;
R_om(2, 2, :) = c_om;
R_om(3, 3, :) = 1;

R_i = zeros(3, 3, length(i));
c_i = cos(i);
s_i = sin(i);
R_i(1, 1, :) = 1;
R_i(2, 2, :) = c_i;
R_i(2, 3, :) = s_i;
R_i(3, 2, :) = -s_i;
R_i(3, 3, :) = c_i;

R_OM = zeros(3, 3, length(OM));
c_OM = cos(OM);
s_OM = sin(OM);
R_OM(1, 1, :) = c_OM;
R_OM(1, 2, :) = s_OM;
R_OM(2, 1, :) = -s_OM;
R_OM(2, 2, :) = c_OM;
R_OM(3, 3, :) = 1;

% Compute combined rotation matrix
R313 = pagemtimes(R_om, R_i); % ECI --> PF
R313 = pagemtimes(R313, R_OM);

% Transpose the combined rotation matrix for PF to ECI transformation
R313_t = pagetranspose(R313);

% Transform position and velocity from Perifocal to ECI
r_ECI = pagemtimes(R313_t, r_PF);
rr = squeeze(r_ECI(:, :, 1));

v_ECI = pagemtimes(R313_t, v_PF);
vv = squeeze(v_ECI(:, :, 1));
end
