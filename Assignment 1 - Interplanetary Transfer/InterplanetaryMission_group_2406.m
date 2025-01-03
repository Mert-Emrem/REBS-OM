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
% First time windows
% dep_window_start   = date2mjd2000([2030, 1, 1, 0, 0, 0]);
% dep_window_end     = date2mjd2000([2044, 4, 1, 0, 0, 0]);
% 
% flyby_window_start = date2mjd2000([2030, 1, 1, 0, 0, 0]);
% flyby_window_end   = date2mjd2000([2044, 4, 1, 0, 0, 0]);
% 
% arr_window_start   = date2mjd2000([2030, 1, 1, 0, 0, 0]);
% arr_window_end     = date2mjd2000([2044, 4, 1, 0, 0, 0]);

%% 2 attempt
% Second time windows

dep_window_start   = date2mjd2000([2040, 1, 1, 0, 0, 0]);
dep_window_end     = date2mjd2000([2044, 4, 1, 0, 0, 0]);

flyby_window_start = date2mjd2000([2040, 1, 1, 0, 0, 0]);
flyby_window_end   = date2mjd2000([2044, 4, 1, 0, 0, 0]);

arr_window_start   = date2mjd2000([2040, 1, 1, 0, 0, 0]);
arr_window_end     = date2mjd2000([2044, 4, 1, 0, 0, 0]);

resolution_dep   = 500; 
resolution_flyby = 500; 
resolution_arr   = 500; 

dep_window   = linspace(dep_window_start,   dep_window_end,   resolution_dep);   % [1 x L]
flyby_window = linspace(flyby_window_start, flyby_window_end, resolution_flyby); % [1 x M]
arr_window   = linspace(arr_window_start,   arr_window_end,   resolution_arr);   % [1 x N]

%% Compute Delta-V Grids for Time Windows
% This step calculates delta-V for transfer arcs in the mission using Lambert's problem
% Create the two matrices containing the DeltaV matrices ( N x M x L )

[M1, M2, Vinf_minus, Vinf_plus] = LambertArcsDeltaV_calculator(dep_window, flyby_window, arr_window,data,...
                                0,... Flag PorkchopPlot
                                0);%  Flag 3Dofs Plot


%% Optimization Using Grid Search + fmincon
% The optimization minimizes delta-V using a hybrid approach combining:
% 1) Grid search to identify approximate solutions
% 2) fmincon for fine-tuning the local minimum


[x_GSf, dv_GSf] = Hybrid_GridSearch_Fmincon(dep_window, flyby_window, arr_window,...
                        M1, M2, Vinf_minus, Vinf_plus, data,...
                        5,... N. trials to run the fmincon
                        0,... Flag Plot
                        0);%  Flag Animated Plot


%% Optimizer: ga
% The optimization minimizes delta-V using ga starting from a refined
% choice of time windows due to the previous results

dep_window_start   = date2mjd2000([2041, 0, 0, 0, 0, 0]);
dep_window_end     = date2mjd2000([2043, 0, 0, 0, 0, 0]);

flyby_window_start = date2mjd2000([2042, 1, 1, 0, 0, 0]);
flyby_window_end   = date2mjd2000([2043, 4, 1, 0, 0, 0]);

arr_window_start   = date2mjd2000([2043, 0, 0, 0, 0, 0]);
arr_window_end     = date2mjd2000([2044, 4, 1, 0, 0, 0]);


resolution_dep   = 700; 
resolution_flyby = 600; 
resolution_arr   = 500;

dep_window   = linspace(dep_window_start,   dep_window_end,   resolution_dep);   % [1 x L]
flyby_window = linspace(flyby_window_start, flyby_window_end, resolution_flyby); % [1 x M]
arr_window   = linspace(arr_window_start,   arr_window_end,   resolution_arr);   % [1 x N]

[x_ga, dv_ga] = GeneticAlgorithm(dep_window, flyby_window, arr_window, data,...
                        20,... N. trials to run the ga
                        1,... Flag Plot
                        0);%  Flag Animated Plot

% fmincon for fine-tuning the minimum
[x_fm, dv_fm] = FminCon(dep_window, flyby_window, arr_window, x_ga, data,...
                        4,... N. trials to run the ga
                        1,... Flag Plot
                        0);%  Flag Animated Plot

%% Optimal Values
% load the optimal values obtained before without running the code
load('Optimal_Values_Struct.mat')

%% Plot for optimal DeltaV: GridSearch + fmincon
% Results with 1000X1000X1000 points (first time windows)
% Same result obtain with 200X200X200 points (second time windows)
% fmincon search with 5 trials running...
%------------------------------------------------------------------
% Optimized departure date from Mercury: 2041 3 29 17 8 2.632399e+00 
% Optimized flyby date via Mars: 2042 6 8 7 47 2.151245e+00 
% Optimized arrival date to Harmonia: 2043 8 5 12 19 2.450843e+01
%-------------------------------------------------------------
% OPTIMAL DELTA v
% Elapsed time is 18.635472 seconds.
%    23.8371
% dep_GS = date2mjd2000([2041 3 29 17 8 2.632399e+00]);
% flyby_GS = date2mjd2000([2042 6 8 7 47 2.151245e+00]);
% arr_GS = date2mjd2000([2043 8 5 12 19 2.450843e+01]);

x = DataOpt.GSf.x;

DeltaV_calculator(x, data, 1);
plotTransfer(x);
Animated_Transfers_Plot(x);


%% Plot for optimal DeltaV: ga
% Optimization finished: average change in the fitness value less than
% options.FunctionTolerance and constraint violation is less than options.ConstraintTolerance.
%---------------------------------------------------------------------
% Optimized departure date from Mercury: 2041 3 30 4 20 3.504044e+01 
% Optimized flyby date via Mars: 2042 6 13 13 43 1.008608e+01 
% Optimized arrival date to Harmonia: 2043 8 21 14 15 3.604346e+01 
%-------------------------------------------------------------
% OPTIMAL DELTA v
% Elapsed time is 1092.091359 seconds.
%    23.8173

x = DataOpt.ga.x;

DeltaV_calculator(x, data, 1)
plotTransfer(x)
Animated_Transfers_Plot(x)






