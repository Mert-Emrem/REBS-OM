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

%%
% Starts and ends of time windows, defined in calendar dates
% (YYYY, M, D, h, m, s)

%% First time windows

% dep_window_start   = date2mjd2000([2030, 1, 1, 0, 0, 0]);
% dep_window_end     = date2mjd2000([2044, 4, 1, 0, 0, 0]);
% 
% flyby_window_start = date2mjd2000([2030, 1, 1, 0, 0, 0]);
% flyby_window_end   = date2mjd2000([2044, 4, 1, 0, 0, 0]);
% 
% arr_window_start   = date2mjd2000([2030, 1, 1, 0, 0, 0]);
% arr_window_end     = date2mjd2000([2044, 4, 1, 0, 0, 0]);

%% Second time windows

dep_window_start   = date2mjd2000([2040, 1, 1, 0, 0, 0]);
dep_window_end     = date2mjd2000([2044, 4, 1, 0, 0, 0]);

flyby_window_start = date2mjd2000([2040, 1, 1, 0, 0, 0]);
flyby_window_end   = date2mjd2000([2044, 4, 1, 0, 0, 0]);

arr_window_start   = date2mjd2000([2040, 1, 1, 0, 0, 0]);
arr_window_end     = date2mjd2000([2044, 4, 1, 0, 0, 0]);

%%

resolution_dep   = 500; 
resolution_flyby = 500; 
resolution_arr   = 500; 

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
                                0,... Flag PorkchopPlot
                                0);%  Flag 3Dofs Plot

%% Optimizer: GridSearch + fmincon
% 1) use the algorithm that sum the two matrices and check if the
% constraint is respected
% 2) Compute the minimum value of DeltaV from the 3D matrix 
% 3) Use this value as first guess of fmnicon to find the local minimum of the
% surface

[x_GSf, dv_GSf] = Hybrid_GridSearch_Fmincon(dep_window, flyby_window, arr_window,...
                        M1, M2, Vinf_minus, Vinf_plus, data,...
                        5,... N. trials to run the fmincon
                        0,... Flag Plot
                        0);%  Flag Animated Plot


%% Plot for optimal DeltaV

% GridSearch + fmincon

% Results with 1000X1000X1000 points (first time windows)
% Same result obtain with 200X200X200 points (second time windows)
% Minimum DeltaV found by the GridSearch Algorithm  ---->23.9193
% fmincon search with 5 trials running...
% RUN number:1
% RUN number:2
% RUN number:3
% RUN number:4
% RUN number:5
% Optimized departure date from Mercury: 2041 3 29 17 8 2.632399e+00 
% Optimized flyby date via Mars: 2042 6 8 7 47 2.151245e+00 
% Optimized arrival date to Harmonia: 2043 8 5 12 19 2.450843e+01 
% Elapsed time is 18.635472 seconds.
%    23.8371
load('x_GSf.mat')
load('dv_GSf.mat')
dep_GS = date2mjd2000([2041 3 29 17 8 2.632399e+00]);
flyby_GS = date2mjd2000([2042 6 8 7 47 2.151245e+00]);
arr_GS = date2mjd2000([2043 8 5 12 19 2.450843e+01]);
plotTransfer([dep_GS, flyby_GS, arr_GS])
DeltaV_calculator([dep_GS, flyby_GS, arr_GS], data, 0)



%% Optimizer: ga

dep_window_start   = date2mjd2000([2041, 0, 0, 0, 0, 0]);
dep_window_end     = date2mjd2000([2043, 0, 0, 0, 0, 0]);

flyby_window_start = date2mjd2000([2042, 1, 1, 0, 0, 0]);
flyby_window_end   = date2mjd2000([2043, 4, 1, 0, 0, 0]);

arr_window_start   = date2mjd2000([2043, 0, 0, 0, 0, 0]);
arr_window_end     = date2mjd2000([2044, 4, 1, 0, 0, 0]);


resolution_dep   = 700; 
resolution_flyby = 600; 
resolution_arr   = 500;

% resolution_dep   = dep_window_end   - dep_window_start;
% resolution_flyby = flyby_window_end - flyby_window_start;
% resolution_arr   = arr_window_end   - arr_window_start;

dep_window   = linspace(dep_window_start,   dep_window_end,   resolution_dep);   % [1 x L]
flyby_window = linspace(flyby_window_start, flyby_window_end, resolution_flyby); % [1 x M]
arr_window   = linspace(arr_window_start,   arr_window_end,   resolution_arr);   % [1 x N]

[x_ga, dv_ga] = GeneticAlgorithm(dep_window, flyby_window, arr_window, data,...
                        20,... N. trials to run the ga
                        1,... Flag Plot
                        0);%  Flag Animated Plot

%% and grid
% Optimization finished: average change in the fitness value less than options.FunctionTolerance and constraint violation is less than options.ConstraintTolerance.
% Optimized departure date from Mercury: 2041 3 30 4 20 3.504044e+01 
% Optimized flyby date via Mars: 2042 6 13 13 43 1.008608e+01 
% Optimized arrival date to Harmonia: 2043 8 21 14 15 3.604346e+01 
% Elapsed time is 1092.091359 seconds.
%    23.8173

load("x_ga.mat")
load("dv_ga.mat")
plotTransfer(x_ga)
% Animated_Transfers_Plot(x_ga)

DeltaV_calculator(x_ga, data, 1)


[x_fm, dv_fm] = FminCon(dep_window, flyby_window, arr_window, x_ga, data,...
                        4,... N. trials to run the ga
                        1,... Flag Plot
                        0);%  Flag Animated Plot

%%