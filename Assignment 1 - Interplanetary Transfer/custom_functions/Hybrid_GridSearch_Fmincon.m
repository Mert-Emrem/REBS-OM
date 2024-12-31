function [x, dv] = Hybrid_GridSearch_Fmincon(dep_window, flyby_window, arr_window, ...
                    M1, M2, Vinf_minus, Vinf_plus, data, NNtrials, FlagPlot, FlagAnimatedPlot)



%% GridSearch Algorithm
% optimizer that combine two deltaV matrices and check if the two lambert
% deltaV values could match together considering the constraint on
% hyperbola physical feasibiliy


Delta_GA = ones(length(flyby_window), length(arr_window), length(dep_window))*NaN;
DeltaVtot = ones(length(flyby_window), length(arr_window), length(dep_window))*NaN;

for i=1:length(flyby_window)
    for j=1:length(arr_window)
        for k=1:length(dep_window)

             if dep_window(k)<flyby_window(i) &&...
                     flyby_window(i)<arr_window(j) &&...
                     not(isnan(M1(i, j, k))) && ...
                     not(isnan(M2(i, j, k)))

             vinf_m = squeeze(Vinf_minus(k,i,:));
             vinf_p = squeeze(Vinf_plus(i,j,:));

             [dvp, ~, rp]  = PowerGravityAssist(vinf_m, vinf_p...
            ,data.Mars.Radius, data.Mars.h_atm, data.Mars.mu, 0);

             Delta_GA(i, j, k) = dvp;

                  if not(isnan(dvp)) && rp>(data.Mars.h_atm+data.Mars.Radius)
                        DeltaVtot(i,j,k) = dvp + M1(i, j, k) + M2(i, j, k);
    
                  end

             end

        end

    end

end


% FIND THE MINUM VALUE OF DELTA_V
[Opt, idx] = min(DeltaVtot(:));
disp(strcat('Minimum DeltaV found by the GridSearch Algorithm  ----> ', num2str(Opt)))
[row, col, depth] = ind2sub(size(DeltaVtot), idx);

% FIND THE BEST 3DOFS OPTIONS FOR fmincon optimizer
t_flyby = flyby_window(row);
t_arr = arr_window(col);
t_dep = dep_window(depth);



%% fmincon


% Save results for each run:
x_trials = zeros([NNtrials, 3]);
dv_trials = zeros([NNtrials, 1]);

disp(['fmincon search with ',num2str(NNtrials),' trials running...']);
tic

x0 = [t_dep, t_flyby, t_arr];

options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...'Display', 'iter', ...
    'MaxIterations', 1000, ...
    'OptimalityTolerance', 1e-10, ...
    'ConstraintTolerance', 1e-6, ...
    'Display','off');

for i = 1:NNtrials

    disp(['RUN number:',num2str(i)]);

    [x, dv] = fmincon(@(x) DeltaV_calculator(x, data, 0), x0, [], [], [], [],...
            [dep_window(1), flyby_window(1), arr_window(1)],...
            [dep_window(end), flyby_window(end), arr_window(end)],...
            @(x) nonlcon(x, data), options);

dv_trials(i, 1) = dv;
x_trials(i, :) = x;

end

depdate = mjd20002date(x(1));
fprintf( ['Optimized departure date from Mercury: ', repmat('%d ', 1, numel(depdate)), '\n'], depdate);

flybydate = mjd20002date(x(2));
fprintf( ['Optimized flyby date via Mars: ', repmat('%d ', 1, numel(depdate)), '\n'], flybydate);

arrdate = mjd20002date(x(3));
fprintf( ['Optimized arrival date to Harmonia: ', repmat('%d ', 1, numel(arrdate)), '\n'], arrdate);


toc

% Select minimum deltaV solution:

[~,index] = min(dv_trials);
dv = dv_trials(index);
disp(dv)
x = x_trials(index,:);

%% plot
if FlagPlot
plotTransfer(x)
end

%% animated plot
if FlagAnimatedPlot
Animated_Transfers_Plot(x)
end

end