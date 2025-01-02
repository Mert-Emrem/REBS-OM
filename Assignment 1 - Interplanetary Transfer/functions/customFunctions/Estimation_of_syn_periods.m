%% Synodic Period and Hohmann Transfer Calculations
%
% DESCRIPTION:
%  This script calculates the synodic periods between celestial bodies 
%  and evaluates transfer times and velocities using Hohmann transfer orbits.
%
% INPUT:
%  - mu_sun: Gravitational parameter of the Sun [km^3/s^2]
%  - a_merc: Semi-major axis of Mercury's orbit [km]
%  - a_mars: Semi-major axis of Mars's orbit [km]
%  - a_harmonia: Semi-major axis of Harmonia's orbit [km]
%
% OUTPUT:
%  - T_syn_1: Synodic period between Mercury and Mars [days]
%  - T_syn_2: Synodic period between Mars and Harmonia [years]
%  - T: Time of flight for Hohmann transfer [s]
%  - T_max: Maximum allowable transfer time [s]
%
% AUTHORS: Richero Giovanni, Emrem Mert, Bernarsconi Ludovico, Serlini
% Mariagiulia
% -------------------------------------------------------------------------

%% Constants
mu_sun = astroConstants(4); % Sun's gravitational parameter [km^3/s^2]

% Semi-major axes of celestial bodies [km]
a_merc = 5.7909e+07;       % Mercury
a_mars = 2.2794e+08;       % Mars
a_harmonia = 3.3922e+08;   % Harmonia (asteroid)

%% Synodic Period Calculations

% Synodic period between Mercury and Mars
T1 = 2 * pi * sqrt(a_merc^3 / mu_sun); % Orbital period of Mercury [s]
T2 = 2 * pi * sqrt(a_mars^3 / mu_sun); % Orbital period of Mars [s]
T_syn_1 = T1 * T2 / abs(T1 - T2);      % Synodic period [s]
T_syn_1 = T_syn_1 / 3600 / 24;         % Convert to days
disp(['Synodic period (Mercury-Mars): ', num2str(T_syn_1), ' days']);

% Synodic period between Mars and Harmonia
T2 = 2 * pi * sqrt(a_mars^3 / mu_sun); % Orbital period of Mars [s]
T3 = 2 * pi * sqrt(a_harmonia^3 / mu_sun); % Orbital period of Harmonia [s]
T_syn_2 = T3 * T2 / abs(T3 - T2);      % Synodic period [s]
T_syn_2 = T_syn_2 / 3600 / 24 / 365;   % Convert to years
disp(['Synodic period (Mars-Harmonia): ', num2str(T_syn_2), ' years']);

%% Hohmann Transfer: Mars to Harmonia

% Velocities of Mars and Harmonia [km/s]
v_mars = sqrt(mu_sun / a_mars);
v_harmonia = sqrt(mu_sun / a_harmonia);

% Semi-major axis of Hohmann transfer orbit [km]
a_h_1 = (a_mars + a_harmonia) / 2;

% Compute time of flight for Hohmann transfer [s]
T = pi * sqrt(a_h_1^3 / mu_sun); % Half of orbital period of transfer ellipse
disp(['Time of flight (Mars-Harmonia): ', num2str(T / 86400), ' days']);

% Maximum allowable transfer time (scaled synodic period) [s]
T_max = 1.2 * T_syn_2 * 365 * 86400; % Convert years to seconds
disp(['Maximum transfer time (Mars-Harmonia): ', num2str(T_max / 86400), ' days']);

%% Hohmann Transfer: Mars to Mercury

% Velocities of Mars and Mercury [km/s]
v_mars = sqrt(mu_sun / a_mars);
v_merc = sqrt(mu_sun / a_merc);

% Semi-major axis of Hohmann transfer orbit [km]
a_h_2 = (a_mars + a_merc) / 2;

% Compute time of flight for Hohmann transfer [s]
T = pi * sqrt(a_h_2^3 / mu_sun); % Half of orbital period of transfer ellipse
disp(['Time of flight (Mars-Mercury): ', num2str(T / 86400), ' days']);
