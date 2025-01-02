function [x, dv] = Hybrid_GridSearch_Fmincon(dep_window, flyby_window, arr_window, ...
                    M1, M2, Vinf_minus, Vinf_plus, data, NNtrials, FlagPlot, FlagAnimatedPlot)
%
% HYBRID_GRIDSEARCH_FMINCON
%
% Function to optimize interplanetary transfer trajectories using a hybrid 
% approach combining Grid Search and fmincon optimization.
%
% INPUT:
%  dep_window        [1xL]    Departure date vector in mjd2000
%  flyby_window      [1xM]    Flyby date vector in mjd2000
%  arr_window        [1xN]    Arrival date vector in mjd2000
%  M1                [MxNxL]  Delta-V matrix for the Mars-Harmonia leg of the transfer
%  M2                [LxMxN]  Delta-V matrix for the Mercury-Mars leg of the transfer
%  Vinf_minus        [LxMx3]  Inbound excess velocity vectors at the flyby [km/s]
%  Vinf_plus         [MxNx3]  Outbound excess velocity vectors at the flyby [km/s]
%  data              [1x1]    Struct including planetary/asteroid parameters 
%                              (radius, atmosphere height, gravitational constant)
%  NNtrials          [int]    Number of initial trials for fmincon optimization
%  FlagPlot          [bool]   Flag to enable plotting of the optimized transfer
%  FlagAnimatedPlot  [bool]   Flag to enable animated visualization of the transfer
%
% OUTPUT:
%  x                 [1x3]    Optimized transfer dates: [departure, flyby, arrival] (mjd2000)
%  dv                [1x1]    Total Delta-V for the optimized trajectory [km/s]
%
% AUTHORS: Richero Giovanni, Emrem Mert, Bernarsconi Ludovico, Serlini
% Mariagiulia
%
% -------------------------------------------------------------------------

% Initialize matrices for Delta-V calculations
Delta_GA = ones(length(flyby_window), length(arr_window), length(dep_window)) * NaN;
DeltaVtot = ones(length(flyby_window), length(arr_window), length(dep_window)) * NaN;

% Grid search loop for Delta-V calculations
for i = 1:length(flyby_window)
    for j = 1:length(arr_window)
        for k = 1:length(dep_window)
            
            % Check date constraints and data validity
            if dep_window(k) < flyby_window(i) && ...
               flyby_window(i) < arr_window(j) && ...
               not(isnan(M1(i, j, k))) && ...
               not(isnan(M2(i, j, k)))

                % Compute Delta-V for the gravity assist at Mars
                vinf_m = squeeze(Vinf_minus(k, i, :));
                vinf_p = squeeze(Vinf_plus(i, j, :));
                [dvp, ~, rp] = PowerGravityAssist(vinf_m, vinf_p, ...
                    data.Mars.Radius, data.Mars.h_atm, data.Mars.mu, 0);

                % Store Delta-V and total Delta-V if constraints are met
                Delta_GA(i, j, k) = dvp;
                if not(isnan(dvp)) && rp > (data.Mars.h_atm + data.Mars.Radius)
                    DeltaVtot(i, j, k) = dvp + M1(i, j, k) + M2(i, j, k);
                end
            end
        end
    end
end

% Identify optimal Delta-V and corresponding indices
[Opt, idx] = min(DeltaVtot(:));
disp(['Minimum DeltaV found by the GridSearch Algorithm ----> ', num2str(Opt)]);
[row, col, depth] = ind2sub(size(DeltaVtot), idx);

% Retrieve optimal dates for fmincon initialization
t_flyby = flyby_window(row);
t_arr = arr_window(col);
t_dep = dep_window(depth);

% Display initial optimal dates
disp('Optimal dates found with Grid Search, first initial guess values for fmincon');
fprintf(['Optimized departure date from Mercury: ', repmat('%d ', 1, numel(mjd20002date(t_dep))), '\n'], mjd20002date(t_dep));
fprintf(['Optimized flyby date via Mars: ', repmat('%d ', 1, numel(mjd20002date(t_flyby))), '\n'], mjd20002date(t_flyby));
fprintf(['Optimized arrival date to Harmonia: ', repmat('%d ', 1, numel(mjd20002date(t_arr))), '\n'], mjd20002date(t_arr));

%% fmincon Optimization

% Initialize storage for multiple trials
x_trials = zeros([NNtrials, 3]);
dv_trials = zeros([NNtrials, 1]);

disp(['fmincon search with ', num2str(NNtrials), ' trials running...']);
tic

x0 = [t_dep, t_flyby, t_arr]; % Initial guess for optimization

% Optimization options
options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'MaxIterations', 1000, ...
    'OptimalityTolerance', 1e-10, ...
    'ConstraintTolerance', 1e-6, ...
    'Display', 'off');

% Run fmincon for each trial
for i = 1:NNtrials
    disp(['RUN number: ', num2str(i)]);
    [x, dv] = fmincon(@(x) DeltaV_calculator(x, data, 0), x0, [], [], [], [], ...
                      [dep_window(1), flyby_window(1), arr_window(1)], ...
                      [dep_window(end), flyby_window(end), arr_window(end)], ...
                      @(x) nonlcon(x, data), options);
    dv_trials(i, 1) = dv;
    x_trials(i, :) = x;
end

% Convert and display optimized dates
depdate = mjd20002date(x(1));
fprintf(['Optimized departure date from Mercury: ', repmat('%d ', 1, numel(depdate)), '\n'], depdate);

flybydate = mjd20002date(x(2));
fprintf(['Optimized flyby date via Mars: ', repmat('%d ', 1, numel(flybydate)), '\n'], flybydate);

arrdate = mjd20002date(x(3));
fprintf(['Optimized arrival date to Harmonia: ', repmat('%d ', 1, numel(arrdate)), '\n'], arrdate);

toc

% Select the solution with minimum Delta-V
[~, index] = min(dv_trials);
dv = dv_trials(index);
disp(['Final optimized Delta-V: ', num2str(dv)]);
x = x_trials(index, :);

%% Plot and Animation
if FlagPlot
    plotTransfer(x);
end
if FlagAnimatedPlot
    Animated_Transfers_Plot(x);
end

end
