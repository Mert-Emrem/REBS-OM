function [T_ode, Y] = Orbit_Analysis(r0, v0, mu, tspan, type_pb, varargin)

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

 



