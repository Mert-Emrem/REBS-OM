function plotTransfer(x)
% Function to plot interplanetary transfer arcs and a Mars flyby
% This function visualizes two transfer arcs given the departure, flyby, 
% and arrival times. Additionally, it plots the flyby trajectory near Mars.
%
% PROTOTYPE:
%  plotTransfer(x)
%
% INPUT:
%  x [3,1]        MJD2000 times array: 
%                 x(1) = departure time (MJD2000)
%                 x(2) = flyby time (MJD2000)
%                 x(3) = arrival time (MJD2000)
%
% OUTPUT:
%  Two visualizations: 
%   1. Interplanetary trajectory with transfer arcs.
%   2. Flyby trajectory around Mars.
%
% AUTHORS: Richero Giovanni, Emrem Mert, Bernarsconi Ludovico, Serlini Mariagiulia
% -------------------------------------------------------------------------

% Extract departure, flyby, and arrival times
t_d = x(1); 
t_f = x(2); 
t_a = x(3); 

%% Compute planetary positions and velocities
[kep_d, ksun] = uplanet(t_d, 1); % Departure (Mercury)
[rr_d, v_merc] = kep2car(kep_d(1), kep_d(2), kep_d(3), kep_d(4), kep_d(5), kep_d(6), ksun);

[kep_f, ksun] = uplanet(t_f, 4); % Flyby (Mars)
[rr_f, v_mars] = kep2car(kep_f(1), kep_f(2), kep_f(3), kep_f(4), kep_f(5), kep_f(6), ksun);

[kep_a, ~, ~] = ephAsteroids(t_a, 40); % Arrival (Harmonia asteroid)
[rr_a, v_harm] = kep2car(kep_a(1), kep_a(2), kep_a(3), kep_a(4), kep_a(5), kep_a(6), ksun);

%% Solve Lambert problems for transfer arcs
[~, ~, ~, ~, vt1_i, VV_minus, ~, ~] = lambertMR(rr_d, rr_f, (t_f - t_d) * 24 * 3600, ksun, 0, 0, 0, 0);
[~, ~, ~, ~, vt2_i, ~, ~, ~] = lambertMR(rr_f, rr_a, (t_a - t_f) * 24 * 3600, ksun, 0, 0, 0, 0);

%% Integrate transfer arcs and planetary orbits
[~, y1] = TwoBodyPb([0 (t_f - t_d) * 24 * 3600], [rr_d; vt1_i'], ksun);
[~, y2] = TwoBodyPb([0 (t_a - t_f) * 24 * 3600], [rr_f; vt2_i'], ksun);

[~, r1] = TwoBodyPb([0 88 * 24 * 3600], [rr_d; v_merc], ksun); % Mercury orbit
[~, r2] = TwoBodyPb([0 690 * 24 * 3600], [rr_f; v_mars], ksun); % Mars orbit
[~, r3] = TwoBodyPb([0 1260 * 24 * 3600], [rr_a; v_harm], ksun); % Harmonia orbit

r1 = r1'; r2 = r2'; r3 = r3';

%% Define celestial objects
Mars.Radius = 3390; Mars.name = 'Mars';
Mercury.Radius = 2439.7; Mercury.name = 'Mercury';
Harmonia.Radius = 107.6; Harmonia.name = 'Harmonia';
Sun.Radius = 700e+3; Sun.name = 'Sun';

%% Plot interplanetary transfers
figure
hold on
plot3(y1(:,1), y1(:,2), y1(:,3), 'LineWidth', 2, 'Color', '#0714fa'); % Lambert arc 1
plot3(y2(:,1), y2(:,2), y2(:,3), 'LineWidth', 2, 'Color', '#e607fa'); % Lambert arc 2
plot3(r1(1,:), r1(2,:), r1(3,:), 'LineStyle','--', 'Color', '#8fe866'); % Mercury orbit
plot3(r2(1,:), r2(2,:), r2(3,:), 'LineStyle','--', 'Color', '#eb8552'); % Mars orbit
plot3(r3(1,:), r3(2,:), r3(3,:), 'LineStyle','--', 'Color', '#bfb588'); % Harmonia orbit

grid on
view(45, 25) % Set viewing angle
axis equal
xlim([-0.41e+9 +0.41e+9]); ylim([-0.41e+9 +0.41e+9]);

% Plot celestial bodies
plotObjects(10, Sun, [0, 0, 0]); % Sun
plotObjects(2000, Mercury, r1(:,1)); % Mercury
plotObjects(1800, Mars, r2(:,1)); % Mars
plotObjects(45000, Harmonia, r3(:,1)); % Harmonia

xlabel('X [Km]', 'fontsize', 14, 'interpreter', 'latex');
ylabel('Y [Km]', 'fontsize', 14, 'interpreter', 'latex');
zlabel('Z [Km]', 'fontsize', 14, 'interpreter', 'latex');
title('Interplanetary Trajectory', 'fontsize', 14, 'interpreter', 'latex');
legend('Lambert arc - 1', 'Lambert arc - 2', 'Mercury - Orbit', 'Mars - Orbit', 'Harmonia - Orbit', 'fontsize', 14, 'interpreter', 'latex');

%% Plot flyby trajectory
vv_inf_minus = VV_minus - v_mars';
vv_inf_plus = vt2_i - v_mars';

mu_mars = astroConstants(14); % Mars gravitational parameter
[~, ~, rp, em, ep, am, ap, vpm, vpp, deltam, deltap] = PowerGravityAssist(vv_inf_minus, vv_inf_plus, Mars.Radius, 220, mu_mars, 1);

YesPlot = ~isnan(rp); % Check for feasible flyby

if YesPlot
    figure
    % time span
    dt = 10; 
    tspan = 0:dt:10000;
    
    rr_p = [rp 0 0]';
    vvpm = [0 vpm 0]';
    y0 = [rr_p; vvpm];
    [~,State] = TwoBodyPb(-tspan,y0,ksun);
    X1 = State(:,1); Y1 = State(:,2); Z1 = State(:,3);
    
    vvpp = [0 vpp 0]';
    y0 = [rr_p; vvpp];
    [~, State] = TwoBodyPb(tspan,y0,ksun);
    X2 = State(:,1); Y2 = State(:,2); Z2 = State(:,3);
    
    
    grid on
    hold on
    view(0, 90) % azimuth and elevation
    axis equal
    xlim([-2e+4 +2e+4]); ylim([-0.41e+5 +0.41e+5]);
    
    X_min_minus = min(X1);
    X_min_plus = min(X2);
    
    
    x0_minus = rp + am;
    x0_plus =  rp  + ap;
    Y_minus = @(x)   tan(pi/2 - deltam /2) * (x-x0_minus);
    Y_plus =  @(x)  -tan(pi/2 - deltap /2)  * (x-x0_plus);

    X = X_min_minus:x0_minus ;
    plot3(X, Y_minus(X), zeros(length(X),1), 'LineWidth', 2, 'Color', '#0714fa')

    X =  X_min_plus:x0_plus ;
    plot3(X, Y_plus(X), zeros(length(X),1), 'LineWidth', 2, 'Color', '#e607fa')

    plot3(X1, Y1, Z1, 'k', 'LineWidth', 2);% Inbound arc 
    plot3(X2, Y2, Z2, 'k', 'LineWidth', 2);% Outbound arc

    plotObjects(1, Mars, [0, 0, 0]); % Mars
    xlabel('X [Km]', 'fontsize', 14, 'interpreter', 'latex');
    ylabel('Y [Km]', 'fontsize', 14, 'interpreter', 'latex');
    title('Flyby at Mars', 'fontsize', 14, 'interpreter', 'latex');
    legend('Inbound arc', 'Outbound arc', 'fontsize', 14, 'interpreter', 'latex');
end
end
