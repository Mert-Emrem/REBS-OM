clc; clear all;


mu_E = astroConstants(13);  % Earth's gravitational parameter [km^3/s^2]
muSun = 1.32712440017987E+11; % Sun's gravitational parameter [km^3/s^2]

resolution = 200;

% Convert departure and arrival dates to MJD2000
date_dep1 = date2mjd2000([2003, 4, 1, 0, 0, 0]);
date_dep2 = date2mjd2000([2003, 8, 1, 0, 0, 0]);
tspan_dep = linspace(date_dep1, date_dep2, resolution);

date_arr1 = date2mjd2000([2003, 9, 1, 0, 0, 0]);
date_arr2 = date2mjd2000([2004, 3, 1, 0, 0, 0]);
tspan_arr = linspace(date_arr1, date_arr2, resolution);

% Initialize the Delta-V matrix
deltaV_totals = NaN(resolution, resolution);

% Precompute ephemeris data for Earth and Mars
earth_state = zeros(resolution, 6);
mars_state = zeros(resolution, 6);

for i = 1:resolution
    % Get ephemerides for Earth and Mars
    ephemeris_dep = uplanet(tspan_dep(i), 3); % Earth
    ephemeris_arr = uplanet(tspan_arr(i), 4); % Mars

    % Convert to Cartesian coordinates
    [r_e, v_e] = kep2car(ephemeris_dep, muSun);
    [r_m, v_m] = kep2car(ephemeris_arr, muSun);

    % Store the state vectors
    earth_state(i, :) = [r_e', v_e'];
    mars_state(i, :) = [r_m', v_m'];
end

% Time of Flight and Lambert's solution
for i = 1:resolution
    for j = 1:resolution
        % Compute the time of flight
        ToF = (tspan_arr(j) - tspan_dep(i)) * 86400; % [s]

        if ToF > 0
            % Solve Lambert's problem
            [~, ~, ~, ERROR, VI, VF, ~, ~] = ...
                lambertMR(earth_state(i, 1:3), mars_state(j, 1:3), ToF, muSun, 0, 0, 0);

            if ERROR == 0
                % Compute Delta-V
                deltaV_1 = vecnorm(VI - earth_state(i, 4:6)); % Departure Delta-V
                deltaV_2 = vecnorm(VF - mars_state(j, 4:6)); % Arrival Delta-V
                deltaV_totals(i, j) = deltaV_1 + deltaV_2;
            end
        end
    end
end

t_offset = datenum('2000-01-01');
[X, Y] = meshgrid(tspan_dep + t_offset, tspan_arr + t_offset); 
figure;
% Subplot 1: 2D Contour Plot
subplot(1, 2, 1);
contourf(X, Y, deltaV_totals', 100, 'LineStyle', 'none'); % Smooth 2D contours
colormap('parula');
clim([5 15]);
colorbar;
xlabel('Departure Date');
ylabel('Arrival Date');
title('2D Porkchop Plot');
datetick('x', 'yyyy-mmm-dd', 'keeplimits'); % Format x-axis as dates
datetick('y', 'yyyy-mmm-dd', 'keeplimits'); % Format y-axis as dates
xtickangle(45);

% Subplot 2: 3D Surface Plot
subplot(1, 2, 2);
surf(X, Y, deltaV_totals', 'EdgeColor', 'none'); % Smooth surface plot
colormap('parula'); % Match colormap with 2D plot
clim([5 15]);
colorbar;
xlabel('Departure Date');
ylabel('Arrival Date');
zlabel('Delta-V (km/s)');
title('3D Porkchop Plot');
datetick('x', 'yyyy-mmm-dd', 'keeplimits'); % Format x-axis as dates

% Doesn't print x axis values for some reason ??

datetick('y', 'yyyy-mmm-dd', 'keeplimits'); % Format y-axis as dates
view(3); % Adjust to a 3D perspective

% Return X and Y to their original MJD values
[X, Y] = meshgrid(tspan_dep, tspan_arr);

% Objective function for Delta-V optimization
deltaV_objective = @(x) interp2(X, Y, deltaV_totals', x(1), x(2), 'spline');

% Initial guess: Use the point corresponding to minimum Delta-V from the matrix
[M, index_of_min] = min(deltaV_totals,[],"all");
[row, col] = ind2sub(size(deltaV_totals), index_of_min);
x0 = [tspan_dep(col), tspan_arr(row)]; % Initial guess (t_dep, t_arr)

% Constraints for fmincon: valid departure and arrival times
lb = [min(tspan_dep), min(tspan_arr)]; % Lower bounds
ub = [max(tspan_dep), max(tspan_arr)]; % Upper bounds
A = [0 -1]; % Ensure arrival time > departure time
b = 0;

% fmincon for constrained optimization
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
[x_opt, fval] = fmincon(deltaV_objective, x0, A, b, [], [], lb, ub, [], options);

% Display refined solution
fprintf('Optimized departure date (MJD2000): %.2f\n', x_opt(1));
fprintf('Optimized arrival date (MJD2000): %.2f\n', x_opt(2));
fprintf('Optimized Delta-V (km/s): %.4f\n', fval);

% fminunc for unconstrained optimization
options = optimoptions('fminunc', 'Display', 'iter', 'Algorithm', 'quasi-newton');
[x_opt_unc, fval_unc] = fminunc(deltaV_objective, x0, options);

% Display refined solution
fprintf('Optimized Delta-V (km/s): %.4f\n', fval_unc);
