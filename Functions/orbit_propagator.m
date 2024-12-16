function [ T, Y ] = orbit_propagator(r0 , v0 , MU , odeXX , options, m , t )
% orbit_propagator propagates an orbit starting from the initial conditions
% v0, r0 for a specified time.
%
% STRUCTURE:
% [T, Y] = orbit_propagator(r0 , v0 , odeXX , options, m, t)
%
% INPUT:
% r0 [3x1]       Initial position [Km]
% v0 [3x1]       Initial velocity [Km]
% MU [1]         Planet's gravitational parameter [km^3/s^2]
% odeXX [1]      Value thath indicates the ODE solver to be used. Legend:
%                   --> 1 ode45
%                   --> 2 ode113
% options [-]    Options used from the numerical solver
% m[1]           Number of time steps to be performed
% t[1]           Time of end for simulation
%
% OUTPUT:
% T[mx1]         Array with the times at each of the m time steps.
% Y[mx6]         Array with the 6 states at each of the m time steps. That is, 
%                row m corresponds to the state at the m-th time step.
%
% CONTRIBUTORS:
% Ludovico Bernaasconi
%
% VERSIONS
% 23-10-2024: First version
%
% -------------------------------------------------------------------------

% options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );


% Definition of the Initial State Vector
y0      = [r0;v0];

% Calculation of the orbital energy and check fi it is a closed orbit
E       = 0.5 * dot(v0,v0) - MU/norm(r0);
if E>= 0
    error("ERROR! The initial conditions represent an open orbit");
end

% Calculation of the Period and time span
tspan   = linspace(0,t,m);

% Perform the integration
switch odeXX
    case 1  %ode45
        [ T, Y ] = ode45( @(t,y) ode_2bp(t,y,MU), tspan, y0, options );
    case 2  %ode113
        [ T, Y ] = ode113( @(t,y) ode_2bp(t,y,MU), tspan, y0, options );
    otherwise
        error ("ERROR! please select a valid ode solver")

end