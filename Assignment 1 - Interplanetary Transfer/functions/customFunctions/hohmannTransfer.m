% Gravitational parameter of the Sun [km^3/s^2]
mu = astroConstants(4);

% Semi-major axes of Mercury and Mars orbits [km]
a1 = 57.909e6; % Mercury's semi-major axis
a2 = 227.956e6; % Mars' semi-major axis

% Compute deltaV and time of flight for Hohmann transfer from Mercury to Mars
[H1, ToF1] = hohmannTransferDeltaV(a1, a2, mu);

% Retrieve orbital elements for the asteroid (ID 40) within the specified arrival window
[kep, f, ~] = ephAsteroids_vec(arr_window, 40);

% Semi-major axis of the asteroid's orbit [km]
a3 = kep(1);

% Compute deltaV and time of flight for Hohmann transfer from Mars to the asteroid
[H2, ToF2] = hohmannTransferDeltaV(a2, a3, mu);

% Total deltaV for the two Hohmann transfers [km/s]
Hohmann = H1 + H2;

% Display the result
disp(['Total deltaV for the Hohmann transfer: ', num2str(Hohmann), ' km/s']);


% Function: Compute deltaV and time of flight for a Hohmann transfer
function [deltaV, ToF] = hohmannTransferDeltaV(a1, a2, mu)
%
% hohmannTransferDeltaV: Computes the total deltaV and time of flight for a
% Hohmann transfer between two circular orbits.
%
% PROTOTYPE:
%   [deltaV, ToF] = hohmannTransferDeltaV(a1, a2, mu)
%
% INPUT:
%   a1 [1]   Semi-major axis of the initial orbit [km]
%   a2 [1]   Semi-major axis of the target orbit [km]
%   mu [1]   Gravitational parameter of the central body [km^3/s^2]
%
% OUTPUT:
%   deltaV [1]   Total deltaV required for the Hohmann transfer [km/s]
%   ToF [1]      Time of flight for the Hohmann transfer [days]
%
% AUTHORS: Richero Giovanni, Emrem Mert, Bernarsconi Ludovico, Serlini
% Mariagiulia
%
% -------------------------------------------------------------------------

    % Orbital velocities for the initial and target circular orbits
    v1 = sqrt(mu / a1); % Velocity in initial orbit [km/s]
    v2 = sqrt(mu / a2); % Velocity in target orbit [km/s]

    % Velocities in the transfer orbit
    v_transfer1 = sqrt(2 * mu * a2 / (a1 * (a1 + a2))); % At the initial orbit [km/s]
    v_transfer2 = sqrt(2 * mu * a1 / (a2 * (a1 + a2))); % At the target orbit [km/s]

    % DeltaV for the transfer burns
    deltaV1 = abs(v_transfer1 - v1); % Burn to enter the transfer orbit [km/s]
    deltaV2 = abs(v2 - v_transfer2); % Burn to enter the target orbit [km/s]

    % Total deltaV for the Hohmann transfer
    deltaV = deltaV1 + deltaV2;

    % Semi-major axis of the transfer orbit
    a_t = (a1 + a2) / 2;

    % Time of flight for the Hohmann transfer [days]
    ToF = pi * sqrt((a_t^3) / mu) / 3600 / 24;

end
