function [r_ECI, v_ECI] = kep2car(a, e, i, OM, om, th, mu)
% Function to convert Keplerian orbital elements to Cartesian state vectors
% This function calculates the position and velocity vectors in the Earth-Centered 
% Inertial (ECI) frame, given the Keplerian orbital elements.
% 
% PROTOTYPE:
%  [r_ECI, v_ECI] = kep2car(a, e, i, OM, om, th, mu)
% 
% INPUT:
%  a    [1x1]  Semi-major axis of the orbit [km]
%  e    [1x1]  Eccentricity of the orbit [-]
%  i    [1x1]  Inclination of the orbit [rad]
%  OM   [1x1]  Right Ascension of the Ascending Node (RAAN) [rad]
%  om   [1x1]  Argument of perigee [rad]
%  th   [1x1]  True anomaly [rad]
%  mu   [1x1]  Gravitational parameter of the central body [km^3/s^2]
% 
% OUTPUT:
%  r_ECI [3x1]  Position vector in the ECI frame [km]
%  v_ECI [3x1]  Velocity vector in the ECI frame [km/s]
% 
% AUTHORS: Richero Giovanni, Emrem Mert, Bernarsconi Ludovico, Serlini Mariagiulia
% -------------------------------------------------------------------------


    % Calculate the semi-latus rectum
    p = a * (1 - e^2);

    % Compute specific angular momentum
    h = sqrt(p * mu);

    % Calculate the orbital radius
    r = p / (1 + e * cos(th));

    % Position and velocity in the Perifocal frame
    r_PF = r * [cos(th), sin(th), 0]';
    v_PF = (mu / h) * [-sin(th), (e + cos(th)), 0]';

    % Define rotation matrices to transform from Perifocal to ECI
    % Rotation matrix for argument of perigee
    R_om = [cos(om),  sin(om),  0;
           -sin(om),  cos(om),  0;
                0,        0,   1];
    
    % Rotation matrix for inclination
    R_i = [1,       0,        0;
           0,  cos(i),  sin(i);
           0, -sin(i),  cos(i)];
    
    % Rotation matrix for RAAN
    R_OM = [cos(OM),  sin(OM),  0;
           -sin(OM),  cos(OM),  0;
                0,        0,   1];
    
    % Combined rotation matrix for ECI to Perifocal (R313 transformation)
    R313 = R_om * R_i * R_OM;

    % Transform position and velocity from Perifocal to ECI
    r_ECI = R313' * r_PF;
    v_ECI = R313' * v_PF;
end


