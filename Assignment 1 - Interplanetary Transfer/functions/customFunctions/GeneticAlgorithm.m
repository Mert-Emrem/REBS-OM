function [x, dv] = GeneticAlgorithm(dep_window, flyby_window, arr_window, data, NNtrials, FlagPlot, FlagAnimatedPlot)
%
% GeneticAlgorithm: Optimize mission trajectory using genetic algorithms
%
% This function uses a genetic algorithm (GA) to optimize the trajectory
% of a space mission by minimizing the total deltaV. Multiple trials are
% executed to ensure convergence and avoid local minima.
%
% PROTOTYPE:
%   [x, dv] = GeneticAlgorithm(dep_window, flyby_window, arr_window, data, NNtrials, FlagPlot, FlagAnimatedPlot)
%
% INPUT:
%   dep_window       [2,1] Departure time window [MJD2000]
%   flyby_window     [2,1] Flyby time window [MJD2000]
%   arr_window       [2,1] Arrival time window [MJD2000]
%   data             Struct containing mission parameters
%   NNtrials         [1]   Number of optimization trials
%   FlagPlot         [1]   Flag to enable transfer plot (1 = plot, 0 = no plot)
%   FlagAnimatedPlot [1]   Flag to enable animated transfer plot (1 = animate, 0 = no animation)
%
% OUTPUT:
%   x                [3,1] Optimized times (departure, flyby, arrival) [MJD2000]
%   dv               [1]   Minimum deltaV of the mission [km/s]
%
% AUTHORS: Richero Giovanni, Emrem Mert, Bernarsconi Ludovico, Serlini
% Mariagiulia
%
% -------------------------------------------------------------------------

% Initialize arrays to store results from multiple trials
x_trials = zeros([NNtrials, 3]); % Optimized times for each trial
dv_trials = zeros([NNtrials, 1]); % DeltaV for each trial

% Display optimization initialization
disp(['GA search with ', num2str(NNtrials), ' trials running...']);
tic

% Number of variables (departure, flyby, arrival times)
nvars = 3; 
% Lower and upper bounds for the variables
lb = [dep_window(1), flyby_window(1), arr_window(1)];  
ub = [dep_window(end), flyby_window(end), arr_window(end)]; 

% GA options
options = optimoptions('ga', ...
    'PopulationSize', 800, ...   % Larger population for better exploration
    'MaxGenerations', 400, ...   % Increase generations for better convergence
    'Display', 'iter');          % Display iteration progress

% Perform multiple GA trials to ensure robust results
for i = 1:NNtrials
    disp(['RUN number: ', num2str(i)]);
    
    % Execute the genetic algorithm
    [x, dv] = ga(@(x) DeltaV_calculator(x, data, 1), nvars, [], [], [], [], ...
                 lb, ub, @(x) nonlcon(x, data), options);

    % Save results for each trial
    dv_trials(i, 1) = dv;
    x_trials(i, :) = x;
end

% Convert optimized times to readable dates and display results
depdate = mjd20002date(x(1));
fprintf(['Optimized departure date from Mercury: ', repmat('%d ', 1, numel(depdate)), '\n'], depdate);

flybydate = mjd20002date(x(2));
fprintf(['Optimized flyby date via Mars: ', repmat('%d ', 1, numel(flybydate)), '\n'], flybydate);

arrdate = mjd20002date(x(3));
fprintf(['Optimized arrival date to Harmonia: ', repmat('%d ', 1, numel(arrdate)), '\n'], arrdate);

toc

% Select the solution with the minimum deltaV across all trials
[~, index] = min(dv_trials);
dv = dv_trials(index);
disp(['Minimum deltaV: ', num2str(dv)]);
x = x_trials(index, :);

%% Plot the optimized transfer trajectory if requested
if FlagPlot
    plotTransfer(x);
end

%% Plot an animated trajectory if requested
if FlagAnimatedPlot
    Animated_Transfers_Plot(x);
end

end
