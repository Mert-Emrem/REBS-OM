function [deltaV_min, ToF] = lambert_deltaV_min(porkchopArg, planet_dep, planet_arr, window_dep, window_arr, deltaV_launcher)

% TODO: Also work in the deltaV_launcher

% Lambert transfer arc with the minimum deltaV given a launch window
% and an optional arrival window
%
% PROTOTYPE:
% [deltaV_min, ToF] = lambert_deltaV_min(planet_dep, planet_arr, ...
%                       window_dep, window_arr, deltaV_max)
%
% INPUT:
% porkchopArg[1] If 1, outputs the porkchop plot, 0 for no plot
% planet_dep[1] Planet index of the departure planet [-]
% planet_arr[1] Planet index of the arrival planet  [-]
%               ...(1 through 8, starting from mercury) 
%
% window_dep[2x6] Departure window (first row: first possible dep.
%                 ...second row: last possible dep.)
%                 [(Gregorian) year, month, day, hour, minute, second]
%
%   **OPTIONAL**
% window_arr[2x6] Arrival window (first row: first possible arr.
%                 ...second row: last possible arr.)
%                 [(Gregorian) year, month, day, hour, minute, second]
%
% deltaV_max[1] Maximum deltaV constraint [km/s]
%
% OUTPUT:
% deltaV_min = The minimum deltaV
% ToF = Time of flight of the transfer that requires the minimum deltaV
%
% Bonus: a plot of the minimum deltaV transfer arc
%
% CONTRIBUTORS:
% Mert Emrem
%
% VERSIONS
% 2024-11-20: 1.0
%

mu_Sun = astroConstants(4); % Sun's gravitational parameter [km^3/s^2]

resolution = 1000;

% Convert departure and arrival dates to MJD2000
date_dep1 = date2mjd2000(window_dep(1,:));
date_dep2 = date2mjd2000(window_dep(2,:));
tspan_dep = linspace(date_dep1, date_dep2, resolution);

if nargin < 5
    % If arrival window not provided, it is set to:
    % from earliest launch to 20 years after last launch
    date_arr1 = date2mjd2000(window_dep(1,:));
    date_arr2 = date2mjd2000(window_dep(2,:)) + 365*20;
else
    date_arr1 = date2mjd2000(window_arr(1,:));
    date_arr2 = date2mjd2000(window_arr(2,:));
end

tspan_arr = linspace(date_arr1, date_arr2, resolution);

% Precompute ephemeris data for given planets
departurePlanet_state = NaN(resolution, 6);
arrivalPlanet_state = NaN(resolution, 6);

for i = 1:resolution
    % Compute ephemerides for departure and arrival planets
    % for both time intervals (departure and arrival windows)

    ephemeris_dep = uplanet(tspan_dep(i), planet_dep); % Earth
    ephemeris_arr = uplanet(tspan_arr(i), planet_arr); % Mars

    % Convert to Cartesian coordinates
    [r_d, v_d] = kep2car(ephemeris_dep, mu_Sun);
    [r_a, v_a] = kep2car(ephemeris_arr, mu_Sun);

    % Store the state vectors
    departurePlanet_state(i, :) = [r_d', v_d'];
    arrivalPlanet_state(i, :) = [r_a', v_a'];
end

% Initialize the Delta-V matrix (as a grid 100 [res x res]),
% such that fminunc can solve for the minimum without being tripped up

% TODO: find another way such that fmincon constraints solve this problem

deltaV_totals = ones(resolution, resolution)*100;

% Time of Flight and Lambert's solution (grid of deltaV)
for i = 1:resolution
    for j = 1:resolution
        % Compute the time of flight
        ToF = (tspan_arr(j) - tspan_dep(i)) * 86400; % [s]
        
        % Ignore temporally impossible solutions (arrival before departure)
        if ToF > 0
            % Solve Lambert's problem
            [~, ~, ~, ~, VI, VF, ~, ~] = ...
                lambertMR(departurePlanet_state(i, 1:3), arrivalPlanet_state(j, 1:3), ...
                ToF, mu_Sun, 0, 0, 0);

            % Compute Delta-V
            deltaV_1 = vecnorm(VI - departurePlanet_state(i, 4:6)); % Departure Delta-V
            deltaV_2 = vecnorm(VF - arrivalPlanet_state(j, 4:6)); % Arrival Delta-V
            deltaV_totals(i, j) = deltaV_1+deltaV_2;
            if deltaV_totals(i, j) >= 100
                deltaV_totals(i, j) = 100;
            end
        end
    end
end

    [X, Y] = meshgrid(tspan_dep, tspan_arr);
    
    % Objective function for Delta-V optimization
    deltaV_objective = @(x) interp2(X, Y, deltaV_totals', x(1), x(2), 'spline');
    
    % Initial guess: Use the point corresponding to minimum Delta-V from the matrix
    [M, index_of_min] = min(deltaV_totals,[],"all");
    [row, col] = ind2sub(size(deltaV_totals), index_of_min);
    x0 = [tspan_dep(col), tspan_arr(row)]; % Initial guess (t_dep, t_arr)

    % fmincon for constrained optimization

    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton');
    [t_opt_unc, fval_unc] = fminunc(deltaV_objective, x0, options);
    
    deltaV_min = fval_unc;
    ToF = t_opt_unc;

    % Display refined solution
    depdate = mjd20002date(t_opt_unc(1));
    fprintf( ['Optimized departure date from planet %d: ', repmat('%d ', 1, numel(depdate)), '\n'], planet_dep, depdate);
    arrdate = mjd20002date(t_opt_unc(2));
    fprintf( ['Optimized arrival date to planet %d: ', repmat('%d ', 1, numel(arrdate)), '\n'], planet_arr, arrdate);
    fprintf('Optimized Delta-V (km/s): %.4f\n', fval_unc);

    if porkchopArg == 1
    
        porkchopPlotter(deltaV_totals,tspan_arr,tspan_dep,fval_unc);
    
    end
end

