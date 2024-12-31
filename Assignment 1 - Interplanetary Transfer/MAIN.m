clear 
close all
clc

addpath 'functions'\timeConversion\time\
addpath 'functions'\
addpath 'custom_functions'\

%% Interplanetary Explorer Mission
% Mercury (dep) -> Mars (flyby) -> Harmonia (arr)
% This assignment involves the design of a mission that begins from Mercury, executes
% a flyby at Marsand concludes in Asteroid No. 40 (henceforth referred to as A40). The
% mission is required to begin no earlier than the beginning of 2030 and conclude before the start
% of 2060. The mission is to begin and end at the same velocities as the departure and arrival
% bodies, thus, orbit insertion and deorbit phases (and their respective ∆V s) are not considered

%% Constant definitions

% Gravitational constants of related bodies
mu_sun = astroConstants(4);
mu_merc = astroConstants(11);
mu_mars = astroConstants(14);

data.Mars.Radius = 3390; % Mars radius [km]
data.Mars.mu = mu_mars;
data.Mars.h_atm = 100;

%% Set of time windows
% hint: comment the time section to see the different results

%% Initial date range (coarse discretization)
% See report for reasoning behind departure windows
date_min = date2mjd2000([2030, 1, 1, 0, 0, 0]);
date_max = date2mjd2000([2044, 1, 1, 0, 0, 0]);
dt = 30;

%% Refined date range (finer discretization)
% choice of dates due from the first attempt
date_min = date2mjd2000([2040, 1, 1, 0, 0, 0]);
date_max = date2mjd2000([2044, 4, 1, 0, 0, 0]);
dt = 10;

%% 

time_window = date_min: dt :date_max;

% Departure window in mjd2000
dep_window = time_window; % [1 x L]

% Flyby window in mjd2000
flyby_window = time_window; % [1 x M]

% Arrival window in mjd2000
arr_window = time_window; % [1 x N]

% Flyby window needs not be same as the departure window,
% but it allows for an intuitive visualization

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








