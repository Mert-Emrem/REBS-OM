function dy = eq_motion_cartesian(t, cart_coord, acc_pert_function_cartesian,  mu) %, J2, Re, day)
% ---------------------------------------------------------------------------------------------
% dy = eq_motion_cartesian(t, cart_coord, acc_pert_function_cartesian,  mu) 
% 
% This function gives as output the system of ODE expressed with the 
% equation of motion in the cartesian reference frame of the primary planet, calling the 
% acc_pert_function_cartesian to calculate the perturbations given by the primary
% planet oblateness and Sun gravitational perturbation as third body
% 
%
% INPUTS
% - ODE variables
%   t                               [1x1]    time variable of the ODE system 
%   cart_coord                      [6x1]    column vector of radius and velocity vectors variables of the
%                                            ODE system: [rx [km]; ry [km]; rz [km]; vx [km/s]; vy [km/s]; vz[km/s]]
% - function for acceleration
%   acc_pert_function_cartesian              @(t, cart_coord) acc_pert_function_cartesiam(t, cart_coord, J2, mu, R_planet, date_mjd2000, tilt)
%   mu                              [1x1]    Primary planet gravitational parameter
%
% OUTPUT
%   dy                              [6x1]    system of 6 ODE in the
%                                            unknown: rx, ry, rz, vx, vy, vz, following the cartesian formulation of
%                                            the problem
%
% Author: Serlini Mariagiulia
% Last update: 23/12/2024
%
% ------------------------------------------------------------------------------------------------

% retrive the unknown radius and the velocity vectors
r = cart_coord(1:3);
v = cart_coord(4:6);

% radius magnitude
rnorm = norm(r);

% acceleration in the cartesian reference frame
a_p = acc_pert_function_cartesian(t, cart_coord);  

% assembly of the ODE system
dy = [ v(1);
       v(2);
       v(3);
       (-mu/rnorm^3)*r(1) + a_p(1);
       (-mu/rnorm^3)*r(2) + a_p(2);
       (-mu/rnorm^3)*r(3) + a_p(3)];

end