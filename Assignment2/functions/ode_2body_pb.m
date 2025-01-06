function dy = ode_2body_pb(~, y, mu, J2, Re)
% ------------------------------------------------
% function dy = ode_2body_pb(~, y, mu, J2, Re)
% 
% this function creates the system of ode to solve the 2 body problem, both
% for the non perturbed and perturbed orbits considering Earth oblateness (J2)
% 
% - NON PERTURBED CASE
%   INPUTS
%   t   [1X1]     Time (can be omitted, as the system is autonomous) [T]
%   y   [6x1]     State of the body ( rx, ry, rz, vx, vy, vz ) [ L, L/T ]
%   mu  [1X1]     Gravitational parameter of the primary [L^3/T^2]
%   
%   OUTPUT
%   dy  [6x1]     System of ODE , derivative of the state [ L/T^2, L/T^3 ]
% 
% - PERTURBED CASE
%   ADDITIONAL INPUTS
%   J2  [1x1]     constant for pertubations
%   Re  [1x1]     Radius of the attractor
%   
%   OUTPUT:
%   dy  [6x1]     System of ODE , derivative of the state [ L/T^2, L/T^3 ]
% 
% Authors: Serlini Mariagiulia, Bernasconi Ludovico, Emrem Mert, Richero
%         Giovanni
% Last update: 23/12/2024
%
%---------------------------------------------------------------------

% retrive the position and velocity vector from the state vector
r = y(1:3);
v = y(4:6);
% magnitude of r
rnorm = norm(r);

% non perturbed case
dy = [ v;
       (-mu/rnorm.^3)*r];

% perturbed case
if nargin>3
a_J2 = ((3/2* J2 * mu * Re^2)/(rnorm^4))*...
       [r(1)/rnorm*(5*(r(3)^2)/(rnorm^2)-1);...
        r(2)/rnorm*(5*(r(3)^2)/(rnorm^2)-1);...
        r(3)/rnorm*(5*(r(3)^2)/(rnorm^2)-3)];

dy = [ v(1);
       v(2);
       v(3);
       (-mu/rnorm^3)*r(1) + a_J2(1);
       (-mu/rnorm^3)*r(2) + a_J2(2);
       (-mu/rnorm^3)*r(3) + a_J2(3)];

end