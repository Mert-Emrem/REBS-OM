function [kep_M] = Orbits_calculator(x, data, flag_atm)
% Orbits_calculator: Calculates orbital elements and maneuvers for a mission
%
% Function to compute the orbital elements and state vectors for departure,
% flyby, and arrival bodies during an interplanetary mission. It solves 
% Lambert's problem for the transfer arcs and returns the resulting 
% orbital elements for all trajectory segments.
%
% PROTOTYPE:
%  [kep_M] = Orbits_calculator(x, data, flag_atm)
%
% INPUT:
%  x [3,1]        Times array:
%                    x(1) = departure time [MJD2000]
%                    x(2) = flyby time [MJD2000]
%                    x(3) = arrival time [MJD2000]
%  data           Struct containing mission parameters:
%                    - data.Mars.Radius: Radius of Mars [km]
%                    - data.Mars.h_atm: Atmospheric constraint altitude [km]
%                    - data.Mars.mu: Gravitational parameter of Mars [km^3/s^2]
%  flag_atm       Atmospheric constraint flag:
%                    0 = No atmospheric constraint during flyby
%                    1 = Atmospheric constraint applied during flyby
%
% OUTPUT:
%  kep_M [4x6]    Orbital elements of the transfer segments:
%                   - kep_M(1,:) = Departure orbit elements
%                   - kep_M(2,:) = Flyby incoming orbit elements
%                   - kep_M(3,:) = Flyby outgoing orbit elements
%                   - kep_M(4,:) = Arrival orbit elements
%
% AUTHORS: Richero Giovanni, Emrem Mert, Bernarsconi Ludovico, 
%          Serlini Mariagiulia
%
% -------------------------------------------------------------------------

% Initialize times
t_d = x(1); % Departure time
t_f = x(2); % Flyby time
t_a = x(3); % Arrival time

% Planet IDs for departure, flyby, and arrival bodies
id_d = 1;  % Mercury
id_f = 4;  % Mars
id_a = 40; % Target Asteroid

% Calculate the orbital elements and state vectors of departure body
[kep_d, ksun] = uplanet(t_d, id_d); 
[rr_d, ~] = kep2car(kep_d(1), kep_d(2), kep_d(3), kep_d(4), kep_d(5), kep_d(6), ksun);

% Calculate the orbital elements and state vectors of the flyby body
[kep_f, ksun] = uplanet(t_f, id_f); 
[rr_f, vv_f] = kep2car(kep_f(1), kep_f(2), kep_f(3), kep_f(4), kep_f(5), kep_f(6), ksun);

% Calculate the orbital elements and state vectors of the arrival body
kep_a = ephAsteroids(t_a, id_a); 
[rr_a, ~] = kep2car(kep_a(1), kep_a(2), kep_a(3), kep_a(4), kep_a(5), kep_a(6), ksun);

% Preallocate orbital elements arrays for the transfer arcs
kep1_i = zeros(1, 6); % Incoming orbit elements for first transfer
kep1_f = zeros(1, 6); % Final orbit elements for first transfer
kep2_i = zeros(1, 6); % Incoming orbit elements for second transfer
kep2_f = zeros(1, 6); % Final orbit elements for second transfer

% Solve Lambert's problem for the first transfer leg
[~, ~, ~, err1, vt1_i, vt1_f, ~, ~] = lambertMR(rr_d, rr_f, (t_f - t_d) * 24 * 3600, ksun, 0, 0, 0, 0);
if err1 == 0
    [kep1_i(1:6)] = car2kep(rr_d, vt1_i);
    [kep1_f(1:6)] = car2kep(rr_f, vt1_f);
end

% Solve Lambert's problem for the second transfer leg
[~, ~, ~, err2, vt2_i, vt2_f, ~, ~] = lambertMR(rr_f, rr_a, (t_a - t_f) * 24 * 3600, ksun, 0, 0, 0, 0);
if err2 == 0
    [kep2_i(1:6)] = car2kep(rr_f, vt2_i);
    [kep2_f(1:6)] = car2kep(rr_a, vt2_f);
end

% Assemble orbital elements matrix
kep_M = [kep1_i; kep1_f; kep2_i; kep2_f];

% Ensure Lambert problem solutions are valid
if (err1 == 0) && (err2 == 0)
    % Compute the hyperbolic excess velocities at the flyby
    vinf_m = vt1_f' - vv_f; % Incoming velocity at Mars
    vinf_p = vt2_i' - vv_f; % Outgoing velocity from Mars

    % Compute deltaV for powered gravity assist at Mars
    if ~flag_atm
        dvp = PowerGravityAssist(vinf_m, vinf_p, data.Mars.Radius, 0, data.Mars.mu, 1);
    else
        dvp = PowerGravityAssist(vinf_m, vinf_p, data.Mars.Radius, data.Mars.h_atm, data.Mars.mu, 1);
    end

    if isnan(dvp) && flag_atm
        warning('Atmospheric constraint prevents feasible hyperbolic transfer.');
    end
end

end
