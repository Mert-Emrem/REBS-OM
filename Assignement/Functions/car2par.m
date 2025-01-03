function [a, e, i, OM, om, th] = car2par(rr, vv, mu)
% ----------------------------------------------------------------
%
% function [a, e, i, OM, om, th] = car2par(rr, vv, mu)
%
% Transformation from cartesian coordinates to Keplerian parameters
%
% STRUCTURE:
% [a, e, i, OM, om, th] = car2par(rr, vv, mu)
%
% INPUTS:
% rr         [3x1]     [km]          position vector                
% vv         [3x1]     [km/s]        velocity vector                
% mu         [1x1]     [km^3/s^2]    gravitational parameter          
%            
% OUTPUTS:            
% a          [1x1]     [km]          semi-major axis                
% e          [1x1]     [-]           eccentircity                   
% i          [1x1]     [rad]         inclination                    
% OM         [1x1]     [rad]         RAAN                           
% om         [1x1]     [rad]         pericenter argument            
% th         [1x1]     [rad]         true argument                  
%
% Author: Serlini Mariagiulia
% Last Update: 30/12/2024
% --------------------------------------------------------------------


% vectors norms
rr_norm = norm(rr);
vv_norm = norm(vv);

% Specific energy and a 
E = vv_norm^2/2 - mu/rr_norm;
a = - mu/2/E;

% h
h = cross(rr, vv);
h_norm = norm(h);

% e
e_vect = cross(vv, h)/mu - rr/rr_norm;
e = norm(e_vect);

% i
k = [0 0 1]';
i = acos(h(3)/h_norm);

% N
N = cross(k, h)/norm(cross(k, h));

% OM
if N(2) >= 0
    OM = acos(N(1));
end
if N(2)<0
    OM = 2*pi - acos(N(1));
end

% om
if e_vect(3) >= 0
    om = acos(dot(N,e_vect)/e);
end
if e_vect(3) < 0
    om = 2*pi - acos(dot(N,e_vect)/e);
end

% th
v_r = dot(vv, rr)/rr_norm;

if v_r >= 0
    th = acos(dot(rr,e_vect)/(rr_norm*e));
    if abs(th - 2*pi) < 1e-6
        th = 0;
    end
end
if v_r < 0
    th = 2*pi - acos(dot(rr,e_vect)/(rr_norm*e));
    if abs(th - 2*pi) < 1e-6
        th = 0;
    end
end


end