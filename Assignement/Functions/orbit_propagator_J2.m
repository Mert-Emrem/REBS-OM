function [ T, Y ] = orbit_propagator_J2(r0 , v0 , MU , odeXX , options, m , n, J2 , Re )
% orbit_propagator propagates an orbit starting from the initial conditions
% v0, r0.
%
% STRUCTURE:
% [T, Y] = orbit_propagator_J2(r0 , v0 , MU , odeXX , options, m , n, J2 , Re )
%
% INPUT:
%  - r0         [3x1]       Initial position [Km]
%  - v0         [3x1]       Initial velocity [Km]
%  - MU         [1]         Planet's gravitational parameter [km^3/s^2]
%  - odeXX      [1]         Value thath indicates the ODE solver to be used. Legend:
%                                   --> 1 ode45
%                                   --> 2 ode113
%  - options    [-]         Options used from the numerical solver
%  - m          [1]         Number of time steps to be performed
%  - n          [1]         Number of periods to be simulated
%  - Re         [1]         Earth's Radius [km]
%  - J2         [1]         Harmonic costant
% 
% OUTPUT:
%  - T          [mx1]       Array with the times at each of the m time steps.
%  - Y          [mx6]       Array with the 6 states at each of the m time steps. That is, 
%                           row m corresponds to the state at the m-th time step.
%
% CONTRIBUTORS:
% Ludovico Bernaasconi
%
% -------------------------------------------------------------------------

% Definition of the Initial State Vector
y0      = [r0;v0];

% Calculation of the orbital energy and check fi it is a closed orbit
E       = 0.5 * dot(v0,v0) - MU/norm(r0);
if E>= 0
    error("ERROR! The initial conditions represent an open orbit");
end

% Calculation of the Period and time span
a       = 1/( 2/norm(r0) - dot(v0,v0)/MU );
T       = 2 * pi * sqrt(a^3 / MU);
tspan   = linspace(0,n*T,m);

% Perform the integration
switch odeXX
    case 1  %ode45
        [ T, Y ] = ode45( @(t,y) ode_2bp_J2(t,y,MU,J2,Re), tspan, y0, options );
    case 2  %ode113
        [ T, Y ] = ode113( @(t,y) ode_2bp_J2(t,y,MU,J2,Re), tspan, y0, options );
    otherwise
        error ("ERROR! please select a valid ode solver")

end