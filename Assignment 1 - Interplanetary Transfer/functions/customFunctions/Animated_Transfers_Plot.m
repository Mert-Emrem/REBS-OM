function Animated_Transfers_Plot(x)
% Function to plot animated interplanetary transfer arcs
% This function visualizes two transfer arcs, given the times of
% departure, flyby, and arrival.
% 
% PROTOTYPE:
%  Animated_Transfers_Plot(x)
% 
% INPUT:
%  x [3,1]        MJD2000 times array: 
%                 x(1) = departure time (MJD2000)
%                 x(2) = flyby time (MJD2000)
%                 x(3) = arrival time (MJD2000)
% 
% OUTPUT:
%  Animated plot of three orbits and the transfer arcs
%
% AUTHORS: Richero Giovanni, Emrem Mert, Bernarsconi Ludovico, Serlini Mariagiulia
% -------------------------------------------------------------------------

% Extract departure, flyby, and arrival times from the input array
t_d = x(1); % Departure time
t_f = x(2); % Flyby time
t_a = x(3); % Arrival time

% Compute the orbital elements and state vectors for the planets and the target asteroid
[kep_d, ksun] = uplanet(t_d, 1); % Orbital elements of Mercury at departure
[rr_d, v_merc] = kep2car(kep_d(1), kep_d(2), kep_d(3), kep_d(4), kep_d(5), kep_d(6), ksun);

[kep_f, ksun] = uplanet(t_f, 4); % Orbital elements of Mars at flyby
[rr_f, v_mars] = kep2car(kep_f(1), kep_f(2), kep_f(3), kep_f(4), kep_f(5), kep_f(6), ksun);

[kep_a, ~, ~] = ephAsteroids(t_a, 40); % Orbital elements of asteroid at arrival
[rr_a, v_harm] = kep2car(kep_a(1), kep_a(2), kep_a(3), kep_a(4), kep_a(5), kep_a(6), ksun);

% Solve Lambert's problem to determine transfer velocities
[~, ~, ~, ~, vt1_i, ~, ~, ~] = lambertMR(rr_d, rr_f, (t_f - t_d) * 24 * 3600, ksun, 0, 0, 0, 0);
[~, ~, ~, ~, vt2_i, ~, ~, ~] = lambertMR(rr_f, rr_a, (t_a - t_f) * 24 * 3600, ksun, 0, 0, 0, 0);

% Integrate the trajectories of the planets and asteroid over time
[~, r1] = TwoBodyPb([0 88*24*3600], [rr_d; v_merc], ksun); % Mercury
[~, r2] = TwoBodyPb([0 690*24*3600], [rr_f; v_mars], ksun); % Mars
[~, r3] = TwoBodyPb([0 1260*24*3600], [rr_a; v_harm], ksun); % Asteroid

% Convert trajectories to column vectors for plotting
r1 = r1';
r2 = r2';
r3 = r3';

% Define physical properties for objects
Mars.Radius = 3390;
Mars.name = 'Mars';

Mercury.Radius = 2439.7;
Mercury.name = 'Mercury';

Harmonia.Radius = 107.6;
Harmonia.name = 'Harmonia';

Sun.Radius = 700e+3;
Sun.name = 'Sun';

% Initialize the animated plot
figure;
hold on;
grid on;
view(45, 25); % Set the viewing angle
axis equal;
xlim([-0.39e+9 +0.41e+9]);
ylim([-0.39e+9 +0.41e+9]);

% Add plot labels and title
xlabel('X [Km]', 'fontsize', 14, 'interpreter', 'latex');
ylabel('Y [Km]', 'fontsize', 14, 'interpreter', 'latex');
zlabel('Z [Km]', 'fontsize', 14, 'interpreter', 'latex');
title('Interplanetary Trajectory', 'fontsize', 14, 'interpreter', 'latex');

% Plot the orbits of Mercury, Mars, and the asteroid
plot3(r1(1,:), r1(2, :), r1(3, :), 'LineStyle', '- -', 'Color', '#8fe866'); % Mercury
plot3(r2(1,:), r2(2, :), r2(3, :), 'LineStyle', '- -', 'Color', '#eb8552'); % Mars
plot3(r3(1,:), r3(2, :), r3(3, :), 'LineStyle', '- -', 'Color', '#bfb588'); % Asteroid

% Add Lambert transfer arcs
line_y1 = animatedline(rr_d(1), rr_d(2), rr_d(3), 'color', '#0714fa', 'LineWidth', 2);
line_y2 = animatedline(rr_f(1), rr_f(2), rr_f(3), 'color', '#e607fa', 'LineWidth', 2);

% Add animated lines for the orbits
line1 = animatedline(rr_d(1), rr_d(2), rr_d(3), 'color', '#8fe866', 'LineWidth', 1.5);
line2 = animatedline(rr_f(1), rr_f(2), rr_f(3), 'color', '#eb8552', 'LineWidth', 1.5);
line3 = animatedline(rr_a(1), rr_a(2), rr_a(3), 'color', '#bfb588', 'LineWidth', 1.5);

% Plot the Sun and initialize planet animations
plotObjects(15, Sun, [0, 0, 0]);
[ObjMerc, Merc_X, Merc_Y, Merc_Z] = Plot_Animated_Objects(2000, Mercury, rr_d);
[ObjMars, Mars_X, Mars_Y, Mars_Z] = Plot_Animated_Objects(2000, Mars, rr_f);
[ObjHarm, Harm_X, Harm_Y, Harm_Z] = Plot_Animated_Objects(50000, Harmonia, rr_a);

% Animate the transfer trajectory over time
ToF1 = (t_f - t_d) * 24 * 3600;
ToF2 = (t_a - t_f) * 24 * 3600;
dt = 250000;
tspan = 0:dt:(ToF1 + ToF2);

% Pause for a moment before starting the animation
pause(3.5);

% Initialize storage variables for trajectory points
r1 = [];
r2 = [];
r3 = [];
y = [];

% Loop through time steps to animate the transfers and object orbits
for jj = 1:length(tspan) - 1
    % Determine the transfer arc (first or second)
    if tspan(jj) <= ToF1
        [~, y_temp] = TwoBodyPb([0 tspan(jj+1)], [rr_d; vt1_i'], ksun);
        y = [y_temp(end, 1:3); y];
        addpoints(line_y1, y(1,1), y(1,2), y(1,3));
    else
        [~, y_temp] = TwoBodyPb([ToF1 tspan(jj+1)], [rr_f; vt2_i'], ksun);
        y = [y_temp(end, 1:3); y];
        addpoints(line_y2, y(1,1), y(1,2), y(1,3));
    end

    % Update object positions for animation
    [~, rtemp1] = TwoBodyPb([0 tspan(jj+1)], [rr_d; v_merc], ksun);
    r1 = [rtemp1(end, 1:3); r1];
    [~, rtemp2] = TwoBodyPb([0 tspan(jj+1)], [rr_f; v_mars], ksun);
    r2 = [rtemp2(end, 1:3); r2];
    [~, rtemp3] = TwoBodyPb([0 tspan(jj+1)], [rr_a; v_harm], ksun);
    r3 = [rtemp3(end, 1:3); r3];

    % Update positions for animated objects
    set(ObjMerc, 'XData', r1(1,1) + Merc_X, 'YData', r1(1,2) + Merc_Y, 'ZData', r1(1,3) + Merc_Z);
    set(ObjMars, 'XData', r2(1,1) + Mars_X, 'YData', r2(1,2) + Mars_Y, 'ZData', r2(1,3) + Mars_Z);
    set(ObjHarm, 'XData', r3(1,1) + Harm_X, 'YData', r3(1,2) + Harm_Y, 'ZData', r3(1,3) + Harm_Z);
    
    addpoints(line1, r1(1,1), r1(1,2), r1(1,3));
    addpoints(line2, r2(1,1), r2(1,2), r2(1,3));
    addpoints(line3, r3(1,1), r3(1,2), r3(1,3));
    
    drawnow;
    pause(0.0001);
end

% Add legend to the plot
legend('Mercury Orbit', 'Mars Orbit', 'Harmonia Orbit', 'Lambert Arc (1)', 'Lambert Arc (2)', ...
       'fontsize', 14, 'interpreter', 'latex');

hold off;
end
