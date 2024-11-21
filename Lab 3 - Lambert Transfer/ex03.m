clc; clear all; close all;

mu_Sun = astroConstants(4); % Sun's gravitational parameter [km^3/s^2]

resolution = 200;

% Convert departure and arrival dates to MJD2000
date_dep1 = date2mjd2000([2003, 4, 1, 0, 0, 0]);
date_dep2 = date2mjd2000([2003, 8, 1, 0, 0, 0]);
tspan_dep = linspace(date_dep1, date_dep2, resolution);

date_arr1 = date2mjd2000([2003, 9, 1, 0, 0, 0]);
date_arr2 = date2mjd2000([2004, 3, 1, 0, 0, 0]);
tspan_arr = linspace(date_arr1, date_arr2, resolution);

% Initialize state matrices for Earth and Mars
earth_state = zeros(resolution, 6);
mars_state = zeros(resolution, 6);

for i = 1:resolution
    % Compute ephemerides for Earth and Mars for both time intervals
    % (departure and arrival windows)
    ephemeris_dep = uplanet(tspan_dep(i), 3); % Earth
    ephemeris_arr = uplanet(tspan_arr(i), 4); % Mars

    % Convert to Cartesian coordinates
    [r_e, v_e] = kep2car(ephemeris_dep, mu_Sun);
    [r_m, v_m] = kep2car(ephemeris_arr, mu_Sun);

    % Store the state vectors
    earth_state(i, :) = [r_e', v_e'];
    mars_state(i, :) = [r_m', v_m'];
end

% Initialize the Delta-V matrix
deltaV_totals = zeros(resolution, resolution);

% Time of Flight and Lambert's solution (grid of deltaV)
for i = 1:resolution
    for j = 1:resolution
        % Compute the time of flight
        ToF = (tspan_arr(j) - tspan_dep(i)) * 86400; % [s]
        
        % Ignore temporally impossible solutions (arrival before departure)
        if ToF > 0
            % Solve Lambert's problem
            [~, ~, ~, ~, VI, VF, ~, ~] = ...
                lambertMR(earth_state(i, 1:3), mars_state(j, 1:3), ToF, mu_Sun, 0, 0, 0);

            % Compute Delta-V
            deltaV_1 = vecnorm(VI - earth_state(i, 4:6)); % Departure Delta-V
            deltaV_2 = vecnorm(VF - mars_state(j, 4:6)); % Arrival Delta-V
            deltaV_totals(i, j) = deltaV_1 + deltaV_2; % Total Delta-V
        end
    end
end

%% Optimized Solution 

[X, Y] = meshgrid(tspan_dep, tspan_arr);

% Objective function for Delta-V optimization
deltaV_objective = @(x) interp2(X, Y, deltaV_totals', x(1), x(2), 'spline');

% Initial guess: Use the point corresponding to minimum Delta-V from the matrix
[~, index_of_min] = min(deltaV_totals,[],"all");
[row, col] = ind2sub(size(deltaV_totals), index_of_min);
x0 = [tspan_dep(col), tspan_arr(row)]; % Initial guess (t_dep, t_arr)

% Constraints for fmincon: valid departure and arrival times
lb = [min(tspan_dep), min(tspan_arr)]; % Lower bounds
ub = [max(tspan_dep), max(tspan_arr)]; % Upper bounds
A = [0 -1]; % Ensure arrival time > departure time
b = 0;

% fmincon for constrained optimization
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
[t_opt, fval] = fmincon(deltaV_objective, x0, A, b, [], [], lb, ub, [], options);

% Display refined solution
depdate = mjd20002date(t_opt(1));
fprintf( ['Optimized departure date: ', repmat('%d ', 1, numel(depdate)), '\n'], depdate);
arrdate = mjd20002date(t_opt(2));
fprintf( ['Optimized arrival date: ', repmat('%d ', 1, numel(arrdate)), '\n'], arrdate);
fprintf('Optimized Delta-V (km/s): %.4f\n', fval);

% fminunc for unconstrained optimization
options = optimoptions('fminunc', 'Display', 'iter', 'Algorithm', 'quasi-newton');
[t_opt_unc, fval_unc] = fminunc(deltaV_objective, x0, options);

% Display refined solution
depdate = mjd20002date(t_opt_unc(1));
fprintf( ['Optimized departure date: ', repmat('%d ', 1, numel(depdate)), '\n'], depdate);
arrdate = mjd20002date(t_opt_unc(2));
fprintf( ['Optimized arrival date: ', repmat('%d ', 1, numel(arrdate)), '\n'], arrdate);
fprintf('Optimized Delta-V (km/s): %.4f\n', fval_unc);
    

%% Orbits Plot

% Earth and mars initial and final state vectors (r_vec, v_vec) in cartesian
% dep: departure, arr: arrival, opt: optimized value

ephemeris_e_dep_opt = uplanet(t_opt(1), 3); % Earth
ephemeris_m_dep_opt = uplanet(t_opt(1), 4); % Mars
ephemeris_m_arr_opt = uplanet(t_opt(2), 4); % Mars
[r_e_dep_opt, v_e_dep_opt] = kep2car(ephemeris_e_dep_opt, mu_Sun);
[r_m_dep_opt, v_m_dep_opt] = kep2car(ephemeris_m_dep_opt, mu_Sun);
[r_m_arr_opt, ~          ] = kep2car(ephemeris_m_arr_opt, mu_Sun);

ToF_opt = (t_opt(2) - t_opt(1))*86400;

[~, ~, ~, ~, VI, ~, ~, ~] = ...
                lambertMR(r_e_dep_opt, r_m_arr_opt, ToF_opt, mu_Sun, 0, 0, 0);

%Lambert Arc Plot
ToF_arc = (t_opt(2)-t_opt(1))*86400; % seconds
tspan = linspace( 0, ToF_arc, 200 );

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
y0_arc = [r_e_dep_opt; VI']';

dy = @(t,y) ode_2body_pb(tspan, y, mu_Sun);
[ T_arc, Y_arc ] = ode113(dy, tspan, y0_arc, options);

%% Full orbits

f1 = figure;

f1.Position = [500 0 1200 500]; % Adjust figure window size

% Orbital periods in days
T_Earth = 365.25; % Earth's orbital period [days]
T_Mars = 686.98;  % Mars's orbital period [days]

% Time spans for full orbits
orbit_time_earth = linspace(0, T_Earth * 86400, 1000); % [s]
orbit_time_mars = linspace(0, T_Mars * 86400, 1000);   % [s]

% Earth full orbit
[r_e_orbit, v_e_orbit] = kep2car(uplanet(t_opt(1), 3), mu_Sun);
y0_e_full = [r_e_orbit; v_e_orbit]';
dy = @(t, y) ode_2body_pb(t, y, mu_Sun);
[~, Y_e_full] = ode113(dy, orbit_time_earth, y0_e_full, options);

% Earth motion over transfer time plot
y0_e = [r_e_dep_opt; v_e_dep_opt]';
dy = @(t,y) ode_2body_pb(t, y, mu_Sun);
[ T_e, Y_e ] = ode113(dy, tspan, y0_e, options);

% Earth motion over departure window plot
y0_e_dep = [earth_state(1, 1:3); earth_state(1, 4:6)]';
dy = @(t,y) ode_2body_pb(tspan_dep, y, mu_Sun);
[ T_e_dep, Y_e_dep ] = ode113(dy, tspan_dep*86400, y0_e_dep, options);

% Mars full orbit
[r_m_orbit, v_m_orbit] = kep2car(uplanet(t_opt(1), 4), mu_Sun);
y0_m_full = [r_m_orbit; v_m_orbit]';
[~, Y_m_full] = ode113(dy, orbit_time_mars, y0_m_full, options);

% Mars motion over transfer time plot
y0_m = [r_m_dep_opt; v_m_dep_opt]';
dy = @(t,y) ode_2body_pb(t, y, mu_Sun);
[ T_m, Y_m ] = ode113(dy, tspan, y0_m, options);

% Earth motion over departure window plot
y0_m_arr = [mars_state(1, 1:3); mars_state(1, 4:6)]';
dy = @(t,y) ode_2body_pb(tspan_arr, y, mu_Sun);
[ T_m_arr, Y_m_arr ] = ode113(dy, tspan_arr*86400, y0_m_arr, options);
hold off

% 1 km = ... AU
km_to_au = 6.6845871226706E-9; 

% Lambert Transfer Arc
p_arc = plot3(Y_arc(:, 1) * km_to_au, Y_arc(:, 2) * km_to_au, Y_arc(:, 3) * km_to_au, ...
    'g', 'LineWidth', 5, 'Color', [0.1569    0.5804    0.2118]);
hold on;

% Mars motion over transfer time plot
plot3(Y_m(:, 1) * km_to_au, Y_m(:, 2) * km_to_au, Y_m(:, 3) * km_to_au, ...
    'r', 'LineWidth', 3);

% Mars motion over arrival window plot
p_m_arrival = plot3(Y_m_arr(:, 1) * km_to_au, Y_m_arr(:, 2) * km_to_au, Y_m_arr(:, 3) * km_to_au, ...
    'r', 'LineWidth', 15);
p_m_arrival.Color(4) = 0.2;

% Earth motion over transfer time plot
plot3(Y_e(:, 1) * km_to_au, Y_e(:, 2) * km_to_au, Y_e(:, 3) * km_to_au, ...
    'b', 'LineWidth', 3);

% Earth motion over departure window plot
p_e_departure = plot3(Y_e_dep(:, 1) * km_to_au, Y_e_dep(:, 2) * km_to_au, Y_e_dep(:, 3) * km_to_au, ...
    'b', 'LineWidth', 15);
p_e_departure.Color(4) = 0.2;

% Earth's Full Orbit
plot3(Y_e_full(:, 1) * km_to_au, Y_e_full(:, 2) * km_to_au, Y_e_full(:, 3) * km_to_au, ...
    'b--', 'LineWidth', 3, 'Color', [0.3490    0.4667    0.8510]);

% Mars's Full Orbit
plot3(Y_m_full(:, 1) * km_to_au, Y_m_full(:, 2) * km_to_au, Y_m_full(:, 3) * km_to_au, ...
    'r--', 'LineWidth', 3, 'Color', [0.8784    0.2196    0.3294]);

% Transfer Positions
plot3(0, 0, 0, 'yo', 'MarkerEdgeColor',[0.8    0.8    0.2863], 'MarkerFaceColor', 'y', 'MarkerSize',25);
plot3(r_e_dep_opt(1) * km_to_au, r_e_dep_opt(2) * km_to_au, r_e_dep_opt(3) * km_to_au, ...
    'bo', 'MarkerFaceColor', 'b', 'MarkerSize',15);
plot3(r_m_arr_opt(1) * km_to_au, r_m_arr_opt(2) * km_to_au, r_m_arr_opt(3) * km_to_au, ...
    'ro', 'MarkerFaceColor', 'r', 'MarkerSize',12);

% Axis Labels and Styling
xlabel('X [AU]');
ylabel('Y [AU]');
zlabel('Z [AU]');
legend('Transfer arc', 'Mars motion during transfer', 'Arrival window', 'Earth motion during transfer', 'Departure window');
grid on;
axis equal;
title('Interplanetary Lambert Transfer');
hold off;

%% Porkchop Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute Time-of-Flight (ToF) lines
ToF_matrix = NaN(resolution, resolution);

for i = 1:resolution
    for j = 1:resolution
        ToF_matrix(i, j) = (tspan_arr(j) - tspan_dep(i)); % ToF in days
    end
end

t_offset = datenum('2000-01-01');
[X, Y] = meshgrid(tspan_dep + t_offset, tspan_arr + t_offset); 
f = figure;

% Subplot 1: 2D Contour Plot
subplot(1, 2, 1);
contourf(X, Y, deltaV_totals', 100, 'LineStyle', 'none'); % Smooth 2D contours
clim([fval_unc fval_unc+10]);
colorbar;
xlabel('Departure Date');
ylabel('Arrival Date');
title('2D Porkchop Plot');

% x-axis values only properly display this way
xticks = linspace(min(tspan_dep + t_offset), max(tspan_dep + t_offset), 6); % Set tick positions
set(gca, 'XTick', xticks); % Apply ticks to x-axis
xticklabels = datestr(xticks, 'yyyy-mmm-dd'); % Generate readable date labels
set(gca, 'XTickLabel', xticklabels); % Set the date labels explicitly
xtickangle(45);

datetick('y', 'yyyy-mmm-dd', 'keeplimits'); % Format y-axis as dates
hold on

% Time-of-Flight lines in days, showing every multiple of ToF_levels
ToF_levels = 60;
ToF_lines_min = round(min(ToF_matrix,[],"all")/ToF_levels)*ToF_levels;
ToF_lines_max = round(max(ToF_matrix,[],"all")/ToF_levels)*ToF_levels;
ToF_lines = ToF_lines_min:ToF_levels:ToF_lines_max; 
contour(X, Y, ToF_matrix', ToF_lines, 'k', 'ShowText', 'on'); 

% Round the minimum Delta V to 4 decimal places, add its point on the 2D plot
Vmin_plot_point = sprintf("   %.4f", fval_unc);
plot(t_opt_unc(1)+t_offset,t_opt_unc(2)+t_offset, 'o', 'MarkerFaceColor', 'y')
text(t_opt_unc(1)+t_offset,t_opt_unc(2)+t_offset,Vmin_plot_point,'Color','yellow','FontSize',14)

subplot(1, 2, 2);
surf(X, Y, deltaV_totals', 'EdgeColor', 'none'); % 3D surface plot
colormap('parula'); 
clim([fval_unc fval_unc+10]);
xlabel('Departure Date');
ylabel('Arrival Date');
zlabel('Delta-V (km/s)');
title('3D Porkchop Plot');

xticks = linspace(min(tspan_dep + t_offset), max(tspan_dep + t_offset), 6); % Set tick positions
set(gca, 'XTick', xticks); % Apply ticks to x-axis
xticklabels = datestr(xticks, 'yyyy-mmm-dd'); % Generate readable date labels
set(gca, 'XTickLabel', xticklabels); % Set the date labels explicitly
xtickangle(-45);

datetick('y', 'yyyy-mmm-dd', 'keeplimits'); % Format y-axis as dates
ytickangle(45);
view(3); % Adjust to a 3D perspective

f.Position = [500 500 1200 500]; % Adjust figure window size