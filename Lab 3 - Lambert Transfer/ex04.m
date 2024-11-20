clc; clear; close all;

%% TODO 
% - Modify the lambert_deltaV_min code such that solutions exist only when
%   (deltaV_2)^2 > (deltaV_total) - (deltaV_1)
%   where deltaV_1 is deltaV_launcher (is this necessary?)
% - Plots do not show what two planets the transfer takes place in

% Required to run this script:
% lambert_deltaV_min
% porkchopPlotter
% uplanet
% lambertOrbitPlotter
% .
% . 
% .

% VENUS -----------------------------------------------------------------
window_dep = [2024, 6, 1, 0, 0, 0; 
              2026, 11, 1, 0, 0, 0];

window_arr = [2024, 12, 1, 0, 0, 0; 
              2027, 6, 1, 0, 0, 0];


deltaV_launcher = 0;

[deltaV_min, ToF] = lambert_deltaV_min(1, 3, 2, window_dep, window_arr, deltaV_launcher);
lambertOrbitPlotter(3, 2, ToF)

% MARS  -----------------------------------------------------------------
window_dep = [2025, 8, 1, 0, 0, 0; 
              2031, 1, 1, 0, 0, 0];

window_arr = [2026, 1, 1, 0, 0, 0; 
              2032, 3, 1, 0, 0, 0];

[deltaV_min, ToF] = lambert_deltaV_min(1, 3, 4, window_dep, window_arr, deltaV_launcher);
lambertOrbitPlotter(3, 4, ToF)