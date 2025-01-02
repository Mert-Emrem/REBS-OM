clear 
close all
clc

addpath 'functions'\timeConversion\time\
addpath 'functions'\
addpath 'custom_functions'\

% INTERPLANETARY EXPLORER MISSION
%
% Mercury (dep) -> Mars (flyby) -> Harmonia (arr)
% This assignment involves the design of a mission that begins from Mercury, executes
% a flyby at Marsand concludes in Asteroid No. 40 (henceforth referred to as A40). The
% mission is required to begin no earlier than the beginning of 2030 and conclude before the start
% of 2060. The mission is to begin and end at the same velocities as the departure and arrival
% bodies, thus, orbit insertion and deorbit phases (and their respective ∆V s) are not considered
%
% N°GROUP: 2406 -REBS Project
%
% AUTHORS: Bernasconi Ludovico, Emrem Mert, Richero Giovanni, Serlini
%          Mariagiulia
%
% -------------------------------------------------------------------------

% Constant definitions

% Gravitational constants of related bodies
mu_sun = astroConstants(4);
mu_merc = astroConstants(11);
mu_mars = astroConstants(14);

data.Mars.Radius = 3390; % Mars radius [km]
data.Mars.mu = mu_mars;
data.Mars.h_atm = 100;

% Starts and ends of time windows, defined in calendar dates
% (YYYY, M, D, h, m, s)

dep_window_start   = date2mjd2000([2030, 1, 1, 0, 0, 0]);
dep_window_end     = date2mjd2000([2044, 4, 1, 0, 0, 0]);

flyby_window_start = date2mjd2000([2030, 1, 1, 0, 0, 0]);
flyby_window_end   = date2mjd2000([2044, 4, 1, 0, 0, 0]);

arr_window_start   = date2mjd2000([2030, 1, 1, 0, 0, 0]);
arr_window_end     = date2mjd2000([2044, 4, 1, 0, 0, 0]);

resolution_dep   = 50; 
resolution_flyby = 50; 
resolution_arr   = 50; 

% Use this instead to set resolution to 1 day

% resolution_dep   = dep_window_end   - dep_window_start
% resolution_flyby = flyby_window_end - flyby_window_start
% resolution_arr   = arr_window_end   - arr_window_start

dep_window   = linspace(dep_window_start,   dep_window_end,   resolution_dep);   % [1 x L]
flyby_window = linspace(flyby_window_start, flyby_window_end, resolution_flyby); % [1 x M]
arr_window   = linspace(arr_window_start,   arr_window_end,   resolution_arr);   % [1 x N]

%% Grid construction
% Create the two matrices containing the DeltaV matrices ( N x M x L )

[M1, M2, Vinf_minus, Vinf_plus] = LambertArcsDeltaV_calculator(dep_window, flyby_window, arr_window,data,...
                                1,... Flag PorkchopPlot
                                1);%  Flag 3Dofs Plot

%% Optimizer: GridSearch + fmincon
% 1) use the algorithm that sum the two matrices and check if the
% constraint is respected
% 2) Compute the minimum value of DeltaV from the 3D matrix 
% 3) Use this value as first guess of fmnicon to find the local minimum of the
% surface

[x_GSf, dv_GSf] = Hybrid_GridSearch_Fmincon(dep_window, flyby_window, arr_window,...
                        M1, M2, Vinf_minus, Vinf_plus, data,...
                        5,... N. trials to run the fmincon
                        1,... Flag Plot
                        0);%  Flag Animated Plot

%% Optimizer: ga

[x_ga, dv_ga] = GeneticAlgorithm(dep_window, flyby_window, arr_window, data,...
                        6,... N. trials to run the ga
                        1,... Flag Plot
                        0);%  Flag Animated Plot








