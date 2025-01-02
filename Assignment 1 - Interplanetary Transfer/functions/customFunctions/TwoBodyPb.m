function [t, y] = TwoBodyPb(tspan, y0, mu)

% Brief explanation:
% This function computes the trajectory of a spacecraft or celestial body 
% in a two-body problem using numerical integration of the equations of motion.

% 
% PROTOTYPE:
%  [t, y] = TwoBodyPb(tspan, y0, mu)
%  
% INPUT:
%  tspan        Vector of times for which the solution is calculated, or 
%               a vector specifying the initial and final times [s].
%  y0   [6,1]   Initial state vector:
%               y0(1:3) = Initial position vector [km]
%               y0(4:6) = Initial velocity vector [km/s]
%  mu   [1]     Gravitational parameter of the primary body [km^3/s^2]
% 
% OUTPUT:
%  t    [x,1]   Vector of times at which the solution is computed [s].
%  y    [x,6]   Solution matrix of the ODE integration:
%               y(:,1:3) = Position vectors [km]
%               y(:,4:6) = Velocity vectors [km/s]
% 
% AUTHORS: Richero Giovanni, Emrem Mert, Bernarsconi Ludovico, Serlini Mariagiulia
% -------------------------------------------------------------------------

    % Define ODE solver options for high precision
    options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
    
    % Perform numerical integration using ode113
    [t, y] = ode113(@TwoBodyEquation, tspan, y0, options, mu);

end
