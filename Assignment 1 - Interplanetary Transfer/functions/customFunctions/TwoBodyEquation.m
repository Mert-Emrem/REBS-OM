function [dy] = TwoBodyEquation(~, y, mu)

% Brief explanation:
% This function defines the equations of motion for the two-body problem,
% providing the time derivatives of the state vector (position and velocity) 
% for use in numerical integration.

% 
% PROTOTYPE:
%  [dy] = TwoBodyEquation(~, y, mu)
%  
% INPUT:
%  y  [6,1]     State vector at a given time:
%               y(1:3) = Position vector [km]
%               y(4:6) = Velocity vector [km/s]
%  mu [1]       Gravitational parameter of the primary body [km^3/s^2]
% 
% OUTPUT:
%  dy [6,1]     Time derivative of the state vector:
%               dy(1:3) = Velocity vector [km/s]
%               dy(4:6) = Acceleration vector [km/s^2]
% 
% AUTHORS: Richero Giovanni, Emrem Mert, Bernarsconi Ludovico, Serlini Mariagiulia
% -------------------------------------------------------------------------

    % Initialize derivative vector
    dy = zeros(6, 1);

    % Compute the norm of the position vector
    r = norm(y(1:3));

    % Velocity components (time derivative of position)
    dy(1:3) = y(4:6);

    % Acceleration components (time derivative of velocity)
    dy(4:6) = -mu / r^3 * y(1:3);
end
