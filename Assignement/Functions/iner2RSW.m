function R = iner2RSW(om, th, i, OM)
% -----------------------------------------------------------
% R = iner2RSW(om, th, i, OM)
%
% This function gives as output the rotational matrix from the cartesian
% inertial reference frame of a planet to the RSW (radial - transversal -
% out-of-plane) of the S/C
%
% INPUTS
% om     [1x1]     [rad]     pericenter anomaly
% th     [1x1]     [rad]     true anomaly
% i      [1x1]     [rad]     inclination
% OM     [1x1]     [rad]     Right Ascension of the Ascending Node
%
% OUTPUT
% R      [3x3]      [-]      Rotational matrix from cartesian to RSW reference frame
%
% Author: Serlini Mariagiulia
% Last Update: 23/12/2024
%
% ----------------------------------------------------------

% rotation around z axis of an angle equal to om+th
R_omth = [cos(om+th), sin(om+th), 0; -sin(om+th), cos(om+th), 0; 0, 0, 1];

% rotation around x of an angle i
R_i = [1, 0, 0; 0, cos(i), sin(i); 0, -sin(i), cos(i)];

% rotation arounx z of an angle OM
R_OM = [cos(OM), sin(OM), 0; -sin(OM), cos(OM), 0; 0, 0, 1];

R = R_omth * R_i * R_OM;

end