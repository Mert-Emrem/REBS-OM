function [rr, vv] = par2car (a, e, i, OM, om, th, mu)
% --------------------------------------------------------------------
% [rr, vv] = par2car (a, e, i, OM, om, th, mu)
% 
% Transformation from Keplerian parameters to cartesian coordinates
%
% INPUTS
% a          [1x1]    [km]                 semi-major axis                
% e          [1x1]    [-]                  eccentircity                   
% i          [1x1]    [rad]                inclination                    
% OM         [1x1]    [rad]                RAAN                           
% om         [1x1]    [rad]                pericenter anomaly             
% th         [1x1]    [rad]                true anomaly                   
% mu         [1x1]    [km^3/s^2]           gravitational parameter        
%    
% --------------------------------------------------------------------
% OUTPUTS
% rr         [3x1]    [km]                 position vector               
% vv         [3x1]    [km/s]               velocity vector               
%
% Authors: Serlini Mariagiulia, Bernasconi Ludovico, Emrem Mert, Richero
%         Giovanni
% Last update: 1/11/2024
%
% --------------------------------------------------------------------


% Rotation matrices to rotate from planetocentric equatorial to perifocal and viceversa
R_OM = [cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];

R_i = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];

R_om = [cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];

R_ECI_PF = R_om*R_i*R_OM;
R_PF_ECI = R_OM'*R_i'*R_om';




% rr e vv in perifocal frame
p = a*(1-e^2);
rr_pf_norm = p/(1+e*cos(th));

rr_pf = rr_pf_norm*[cos(th) sin(th) 0]';
vv_pf = sqrt(mu/p)*[-sin(th) e+cos(th) 0]';

% rr e vv in the planetocentric equatorial rf
rr = R_PF_ECI*rr_pf;
vv = R_PF_ECI*vv_pf;

end
