function [x, dv] = GeneticAlgorithm(dep_window, flyby_window, arr_window, data, NNtrials,  FlagPlot, FlagAnimatedPlot)

% To check algorithm convergence set NNtrials > 1.

% Save results for each run:
x_trials = zeros([NNtrials, 3]);
dv_trials = zeros([NNtrials, 1]);

% Display message
disp(['ga search with ',num2str(NNtrials),' trials running...']);
tic

nvars = 3; % Number of variables (departure, flyby, arrival times)
lb = [dep_window(1), flyby_window(1), arr_window(1)];  % Lower bounds
ub = [dep_window(end), flyby_window(end), arr_window(end)]; % Upper bounds

options = optimoptions('ga', ...
    'PopulationSize', 200, ...  % Larger population for better exploration
    'MaxGenerations', 500, ...  % Increase generations for improved convergence
     'Display', 'iter');

for i = 1:NNtrials
    
     disp(['RUN number:',num2str(i)]);
    % Run genetic algorithm
    [x, dv] = ga(@(x) DeltaV_calculator(x, data, 1), nvars, [], [], [], [], ...
                 lb, ub, @(x) nonlcon(x, data), options);

    % Save results
    dv_trials(i, 1) = dv;
    x_trials(i, :) = x;

end



% Convert results to readable dates
depdate = mjd20002date(x(1));
fprintf(['Optimized departure date from Mercury: ', repmat('%d ', 1, numel(depdate)), '\n'], depdate);

flybydate = mjd20002date(x(2));
fprintf(['Optimized flyby date via Mars: ', repmat('%d ', 1, numel(flybydate)), '\n'], flybydate);

arrdate = mjd20002date(x(3));
fprintf(['Optimized arrival date to Harmonia: ', repmat('%d ', 1, numel(arrdate)), '\n'], arrdate);

toc

% Select minimum deltaV solution:
[~, index] = min(dv_trials);
dv = dv_trials(index);
disp(dv)
x = x_trials(index,:);

%% plot
if FlagPlot
% Plot optimal transfer
plotTransfer([x(1), x(2), x(3)])
end

%% animated plot
if FlagAnimatedPlot
Animated_Transfers_Plot(x)
end

end


% options = optimoptions('ga', ...
%     'PopulationSize', 200, ...  % Larger population for better exploration
%     'MaxGenerations', 500, ...  % Increase generations for improved convergence
%     'Display', 'iter', ...
%     'EliteCount', ceil(0.05 * 200), ... % 5% of population as elite
%     'CrossoverFraction', 0.85, ... % Higher crossover to encourage exploitation
%     'MutationFcn', {@mutationadaptfeasible, 0.1}, ... % Higher mutation rate for diversity
%     'PlotFcn', @gaplotbestf, ...
%     'NonlinConAlgorithm','penalty',...
%     'ConstraintTolerance',1e-6,...
%     'FunctionTolerance',1e-8,...
%     'MaxStallGenerations', 30)  % Allow more stalls before termination

    % GA options
% options = optimoptions('ga', ...
%     'PopulationSize', 300, ...               % Larger population for complex constraints
%     'MaxGenerations', 700, ...               % More generations to allow convergence
%     'EliteCount', ceil(0.05 * 300), ...      % 5% elite preservation
%     'CrossoverFraction', 0.85, ...           % High crossover to accelerate convergence
%     'MutationFcn', {@mutationadaptfeasible, 0.12}, ... % Adaptive feasible mutation
%     'FunctionTolerance', 1e-10, ...          % High precision
%     'ConstraintTolerance', 1e-6, ...         % Tight constraint handling
%     'NonlinConAlgorithm', 'auglag', ...      % Augmented Lagrangian for nonlinear constraints
%     'MaxStallGenerations', 40, ...           % Avoid premature convergence
%     'Display', 'iter');
