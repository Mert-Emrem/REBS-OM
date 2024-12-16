function [rr, vv] = par2car (a, e, i, OM, om, th, mu)

% Transformation from Keplerian parameters to cartesian coordinates
%
% [rr, vv] = par2car (a, e, i, OM, om, th, mu)
%
% --------------------------------------------------------------------
% Input argument
% a          [1x1] semi-major axis                [km]
% e          [1x1] eccentircity                   [-]
% i          [1x1] inclination                    [rad]
% OM         [1x1] RAAN                           [rad]
% om         [1x1] pericenter anomaly             [rad]
% th         [1x1] true anomaly                   [rad]
% mu         [1x1] gravitational parameter        [km^3/s^2]   per i progetti mu = 398600 km^3/s^2
%
% --------------------------------------------------------------------
% Output argument
% rr         [3x1] position vector                [km]
% vv         [3x1] velocity vector                [km/s]
%
% --------------------------------------------------------------------


% DEFINIZNIONE MATRICI DI ROTAZIONE
R_OM = [cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];

R_i = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];

R_om = [cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];

R_ECI_PF = R_om*R_i*R_OM;
R_PF_ECI = R_OM'*R_i'*R_om';

% calcolo norme


% calcolo di rr e vv nel sistema perifocale
p = a*(1-e^2);
rr_pf_norm = p/(1+e*cos(th));

rr_pf = rr_pf_norm*[cos(th) sin(th) 0]';
vv_pf = sqrt(mu/p)*[-sin(th) e+cos(th) 0]';

% calcolo vettori rr e vv nel sistema equatoriale geocentrico
rr = R_PF_ECI*rr_pf;
vv = R_PF_ECI*vv_pf;

end
