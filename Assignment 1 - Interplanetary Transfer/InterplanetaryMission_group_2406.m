%% INTERPLANETARY EXPLORER MISSION
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% This script designs an interplanetary mission that begins from Mercury, 
% executes a flyby at Mars, and concludes at Asteroid 40 Harmonia. The mission 
% avoids orbit insertion/deorbit phases by assuming departure and arrival 
% velocities match those of the respective celestial bodies.
%
%--------------------------------------------------------------------------
%
% OBJECTIVES:
% - Optimize the mission trajectory to minimize delta-V
% - Explore multiple time windows using Grid Search and Genetic Algorithms
% - Visualize the optimized trajectories
%
% N°GROUP: 2406 -REBS Project
%
% AUTHORS: Bernasconi Ludovico, Emrem Mert, Richero Giovanni, Serlini
%          Mariagiulia
%
% -------------------------------------------------------------------------

%% Initialization
clear; close all; clc;

% Graphical Orbital Validation 
% All the calculations described were thoroughly checked using Animated_Transfers_Plot.m,
% that takes the 3 DoFs as input’s values and produces an animated
% visualization of the orbits and planet positions with the transfer trajectory, 
% facilitating the identification of optimal maneuver windows. Additionally, the function 
% synchronizes the mission timeline with the time steps used by Matlab’s ode solvers. Since only the
% integration time boundaries, not the individual time steps, can be specified, the function uses a loop
% to compute position and velocity iteratively. This approach guarantees that the transfer trajectory
%is displayed accurately and consistently throughout the simulation.


% Add required paths for auxiliary functions
addpath 'functions'\timeConversion\time\
addpath 'functions'\
addpath 'functions'\'customFunctions'\

% Constant definitions

% Gravitational constants of related bodies
mu_sun = astroConstants(4);
mu_merc = astroConstants(11);
mu_mars = astroConstants(14);
% Struct containing data of Mars
data.Mars.Radius = 3390; % Mars radius [km]
data.Mars.mu = mu_mars;
data.Mars.h_atm = 220;

%% %% Define Mission Time Windows
% Starts and ends of time windows, defined in calendar dates
% (YYYY, M, D, h, m, s)

%% 1 attempt (used but commented)
% % First time windows chosing resolution time window = 100
% dep_window_start   = date2mjd2000([2030, 1, 1, 0, 0, 0]);
% dep_window_end     = date2mjd2000([2044, 4, 1, 0, 0, 0]);
% 
% flyby_window_start = date2mjd2000([2030, 1, 1, 0, 0, 0]);
% flyby_window_end   = date2mjd2000([2044, 4, 1, 0, 0, 0]);
% 
% arr_window_start   = date2mjd2000([2030, 1, 1, 0, 0, 0]);
% arr_window_end     = date2mjd2000([2044, 4, 1, 0, 0, 0]);
%

%% 2 attempt
% Second time windows
dep_window_start   = date2mjd2000([2040, 1, 1, 0, 0, 0]);
dep_window_end     = date2mjd2000([2044, 4, 1, 0, 0, 0]);

flyby_window_start = date2mjd2000([2040, 1, 1, 0, 0, 0]);
flyby_window_end   = date2mjd2000([2044, 4, 1, 0, 0, 0]);

arr_window_start   = date2mjd2000([2040, 1, 1, 0, 0, 0]);
arr_window_end     = date2mjd2000([2044, 4, 1, 0, 0, 0]);

%% 3 attempt
% Third time windows
% dep_window_start   = date2mjd2000([2040, 1, 1, 0, 0, 0]);
% dep_window_end     = date2mjd2000([2042, 1, 1, 0, 0, 0]);
% 
% flyby_window_start = date2mjd2000([2042, 1, 1, 0, 0, 0]);
% flyby_window_end   = date2mjd2000([2043, 1, 1, 0, 0, 0]);
% 
% arr_window_start   = date2mjd2000([2043, 1, 1, 0, 0, 0]);
% arr_window_end     = date2mjd2000([2044, 4, 1, 0, 0, 0]);

%% Time windows
resolution_dep   = 150; 
resolution_flyby = 150; 
resolution_arr   = 150; 

dep_window   = linspace(dep_window_start,   dep_window_end,   resolution_dep);   % [1 x L]
flyby_window = linspace(flyby_window_start, flyby_window_end, resolution_flyby); % [1 x M]
arr_window   = linspace(arr_window_start,   arr_window_end,   resolution_arr);   % [1 x N]

%% Compute Delta-V Grids for Time Windows
% This step calculates delta-V for transfer arcs in the mission using Lambert's problem
% Create the two matrices containing the DeltaV matrices ( N x M x L )

[M1, M2, Vinf_minus, Vinf_plus] = LambertArcsDeltaV_calculator(dep_window, flyby_window, arr_window,data,...
                                1,... Flag PorkchopPlot
                                1);%  Flag 3Dofs Plot


%% Optimization Using Grid Search + fmincon
% The optimization minimizes delta-V using a hybrid approach combining:
% 1) Grid search to identify approximate solutions
% 2) fmincon for fine-tuning the local minimum


[x_GSf, dv_GSf] = Hybrid_GridSearch_Fmincon(dep_window, flyby_window, arr_window,...
                        M1, M2, Vinf_minus, Vinf_plus, data,...
                        5,... N. trials to run the fmincon
                        1,... Flag Plot
                        1);%  Flag Animated Plot


%% Optimizer: ga
% The optimization minimizes delta-V using ga starting from a refined
% choice of time windows due to the previous results

dep_window_start   = date2mjd2000([2041, 0, 0, 0, 0, 0]);
dep_window_end     = date2mjd2000([2043, 0, 0, 0, 0, 0]);

flyby_window_start = date2mjd2000([2042, 1, 1, 0, 0, 0]);
flyby_window_end   = date2mjd2000([2043, 4, 1, 0, 0, 0]);

arr_window_start   = date2mjd2000([2043, 0, 0, 0, 0, 0]);
arr_window_end     = date2mjd2000([2044, 4, 1, 0, 0, 0]);


dep_window   = [dep_window_start,   dep_window_end];
flyby_window = [flyby_window_start, flyby_window_end];
arr_window   = [arr_window_start,   arr_window_end];   % [1 x N]

[x_ga, dv_ga] = GeneticAlgorithm(dep_window, flyby_window, arr_window, data,...
                        10,... N. trials to run the ga
                        0,... Flag Plot
                        0);%  Flag Animated Plot

% fmincon for fine-tuning the minimum
[x_fm, dv_fm] = FminCon(dep_window, flyby_window, arr_window, x_ga, data,...
                        4,... N. trials to run the ga
                        0,... Flag Plot
                        0);%  Flag Animated Plot



x_min = x_fm;
dv_min = dv_fm;

fprintf('Minimum total DeltaV considering all optimization methodologies')
depdate = mjd20002date(x_min(1));
fprintf(['Optimized departure date from Mercury: ', repmat('%d ', 1, numel(depdate)), '\n'], depdate);

flybydate = mjd20002date(x_min(2));
fprintf(['Optimized flyby date via Mars: ', repmat('%d ', 1, numel(depdate)), '\n'], flybydate);

arrdate = mjd20002date(x_min(3));
fprintf(['Optimized arrival date to Harmonia: ', repmat('%d ', 1, numel(arrdate)), '\n'], arrdate);



%% Optimal Values
% load the optimal values obtained before without running the code
% Optimized departure date from Mercury: 2041 3 29 11 27 1.429861e+01 
% Optimized flyby date via Mars: 2042 6 7 0 32 6.662251e+00 
% Optimized arrival date to Harmonia: 2043 7 31 22 44 5.719183e+01 
% Elapsed time is 0.467924 seconds.
%    23.8269

%% Plot for optimal DeltaV: GridSearch + fmincon
% Results with 1000X1000X1000 points (first time windows)
% Same result obtain with 200X200X200 points (second time windows)

DeltaV_calculator(x, data, 1);
plotTransfer(x_min);
Animated_Transfers_Plot(x_min);


%% Optimal DeltaV: ga + fmincon
% Outcome and results:

td = [2041 3 29 11 27 1.429861e+01];
tf = [2042 6 7   0 32 6.662251e+00 ];
ta = [2043 7 31 22 44 5.719183e+01];
td = date2mjd2000(td);
tf = date2mjd2000(tf);
ta = date2mjd2000(ta);

x = [td, tf, ta];

%%
DeltaV_calculator(x, data, 1)
%%
plotTransfer(x)
%%
M = Orbits_calculator(x, data, 1);
%%
Animated_Transfers_Plot(x)









