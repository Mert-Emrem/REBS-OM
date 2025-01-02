function [x, dv] = FminCon(dep_window, flyby_window, arr_window, x0, data, NNtrials, FlagPlot, FlagAnimatedPlot)
%
% FminCon: Optimize mission trajectory using fmincon
%
% This function performs trajectory optimization using `fmincon` to minimize
% the total deltaV of the mission. Multiple trials are run with random
% initial guesses to avoid local minima.
%
% PROTOTYPE:
%   [x, dv] = FminCon(dep_window, flyby_window, arr_window, x0, data, NNtrials, FlagPlot, FlagAnimatedPlot)
%
% INPUT:
%   dep_window     [2,1] Departure time window [MJD2000]
%   flyby_window   [2,1] Flyby time window [MJD2000]
%   arr_window     [2,1] Arrival time window [MJD2000]
%   x0             [3,1] Initial guess for times (departure, flyby, arrival)
%   data           Struct containing mission parameters
%   NNtrials       [1]   Number of optimization trials
%   FlagPlot       [1]   Flag to enable transfer plot (1 = plot, 0 = no plot)
%   FlagAnimatedPlot [1] Flag to enable animated transfer plot (1 = animate, 0 = no animation)
%
% OUTPUT:
%   x              [3,1] Optimized times (departure, flyby, arrival) [MJD2000]
%   dv             [1]   Minimum deltaV of the mission [km/s]
%
% AUTHORS: Richero Giovanni, Emrem Mert, Bernarsconi Ludovico, Serlini
% Mariagiulia
%
% -------------------------------------------------------------------------

% Initialize arrays to store results of multiple trials
x_trials = zeros([NNtrials, 3]); % Optimized times for each trial
dv_trials = zeros([NNtrials, 1]); % DeltaV for each trial

% Display optimization progress
disp(['fmincon search with ', num2str(NNtrials), ' trials running...']);
tic

% Optimization options
options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ... % Sequential Quadratic Programming algorithm
    'MaxIterations', 1000, ...
    'OptimalityTolerance', 1e-10, ...
    'ConstraintTolerance', 1e-6, ...
    'Display', 'off'); % Suppress output

% Loop through multiple trials to avoid local minima
for i = 1:NNtrials
    disp(['RUN number: ', num2str(i)]);

    % Run fmincon optimization
    [x, dv] = fmincon(@(x) DeltaV_calculator(x, data, 0), ... % Objective function
                      x0, [], [], [], [], ... % Initial guess and no linear constraints
                      [dep_window(1), flyby_window(1), arr_window(1)], ... % Lower bounds
                      [dep_window(end), flyby_window(end), arr_window(end)], ... % Upper bounds
                      @(x) nonlcon(x, data), ... % Nonlinear constraints
                      options);

    % Store results
    dv_trials(i, 1) = dv;
    x_trials(i, :) = x;
end

% Convert optimized times to readable dates and display
depdate = mjd20002date(x(1));
fprintf(['Optimized departure date from Mercury: ', repmat('%d ', 1, numel(depdate)), '\n'], depdate);

flybydate = mjd20002date(x(2));
fprintf(['Optimized flyby date via Mars: ', repmat('%d ', 1, numel(depdate)), '\n'], flybydate);

arrdate = mjd20002date(x(3));
fprintf(['Optimized arrival date to Harmonia: ', repmat('%d ', 1, numel(arrdate)), '\n'], arrdate);

toc

% Select the solution with the minimum deltaV from all trials
[~, index] = min(dv_trials);
dv = dv_trials(index);
disp(dv)
x = x_trials(index, :);

%% Plot transfer trajectory if requested
if FlagPlot
    plotTransfer(x)
end

%% Plot animated trajectory if requested
if FlagAnimatedPlot
    Animated_Transfers_Plot(x)
end

end
