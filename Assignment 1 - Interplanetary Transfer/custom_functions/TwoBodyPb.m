function [t,y] = TwoBodyPb(tspan,y0,mu)
% TwoBodyPb: Two-body problem integration 
% 
% Function to compute trajectory in the two-body problem.
% 
% PROTOTYPE:
%  [t,y] = 2BodyPb(tspan,y0,mu)
%  
% INPUT:
%  tspan        Vector of times in which the solution is calculated or 
%               vector of initial and final times [s]
%  y0  [6,1]    Initial condition vector [km] and [km/s]
%  mu  [1]      Gravitational parameter of the primary body [km^3/s^2]
% 
% OUTPUT:
%  t  [x,1]     Vector of times in which the solution is calculated [s]
%  y  [x,6]     Solution matrix of the ode integration [km] and [km/s]

options = odeset('RelTol',1e-13,'AbsTol',1e-14);
[t,y] = ode113(@TwoBodyEquation,tspan,y0,options,mu);

end