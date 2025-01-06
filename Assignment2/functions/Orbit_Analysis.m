function [T_ode, Y] = Orbit_Analysis(r0, v0, mu, tspan, type_pb, varargin)
% ------------------------------------------------------------------
% function [T_ode, Y] = Orbit_Analysis(r0, v0, mu, tspan, type_pb, varargin)
%
% This function propagates an unperturbed or perturbed orbit (only Earth
% oblatenes perturbation) calling the ode_2body_pb to create the ODE system
% and using ode113 to solve it
%
% INPUTS
% r0        [3x1]         [km]             initial position vector as column vector
% v0        [3x1]         [km/s]           initial velocity vector as column vector
% mu        [1x1]         [km^3/s^2]       gravitational constant of the primary planet
% tsapn     [nx1]         [s]              time window in which to solve the problem
% type_pb   [chart]       [-]              specify the unperturbed or pertubed problem
%                                          - 'non_perturbed'
%                                          - 'perturbed'
%     
% ADDITIONAL INPUTS if 'perturbed'     
% J2        [1x1]         [-]              Gravitational Harmonic coefficent of the primary planet
% R         [1x1]         [km]             radius of the primary planet
%       
% OUTPUTS       
% T_ode     [nx1]         [s]              vector of time in which the ODE system have been solved
% Y         [nx1]         [km, km/s]       Solution f the ode problem
%                                          col 1-3: rx, ry, rz
%                                          col 4-6: vx, vy, vz
%     
% Authors: Serlini Mariagiulia, Bernasconi Ludovico, Emrem Mert, Richero
%         Giovanni
% Last update: 23/12/2024
%
% -------------------------------------------------------------------

% state vector
y0 = [r0; v0]; 

% Set options for the ODE solver
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 

if strcmp(type_pb, 'non_perturbed')
        [T_ode, Y] = ode113(@(t,y) ode_2body_pb(t, y, mu), tspan, y0, options);  
    elseif strcmp(type_pb, 'perturbed')
        J2 = varargin{1};
        Re = varargin{2};
        [T_ode, Y] = ode113(@(t,y) ode_2body_pb(t, y, mu, J2, Re), tspan, y0, options);    
    else
        error('Not recognised. Use ''non_perturbed'' or ''perturbed''.')
end

 



