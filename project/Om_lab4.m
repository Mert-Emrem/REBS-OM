clear all
close all
clc

%% Exercise 1a

mu_Sun = astroConstants(4);
mu_Earth = astroConstants(13);
AU = astroConstants(2);

v_inf_minus = [15.1; 0; 0]; % VF Merc-Mars
v_inf = norm(v_inf_minus);
Delta = 9200;
r = [1; 0; 0]*AU;

V_Earth = sqrt(mu_Sun/AU);
V_Earth= V_Earth * [0 1 0]';

a = - mu_Earth/v_inf.^2;
turnAng = 2*atan2(-a, Delta);
e = 1/sin(turnAng/2);
r_p = a*(1-e);
delta_V = 2*v_inf*sin(turnAng/2);

VV_minus = V_Earth+ v_inf_minus;


assist_types = {'leading', 'trailing', 'under', 'over'};
assist_type = assist_types{2};



norm2helio = cross(r, VV_minus);
norm2helio = norm2helio/norm(norm2helio);

    switch assist_type
        case 'leading'
            rot_dir = -norm2helio;
        case 'trailing'
            rot_dir = +norm2helio;
        case 'under'
            rot_dir = -V_Earth;      
    end
rot_dir = rot_dir/norm(rot_dir);

v_inf_plus = RodrigueRotationVector(v_inf_minus, turnAng, rot_dir);
VV_plus = V_Earth+ v_inf_plus;




figure

% plot of the sun (too small compared to the heliptic arc)
% ... I have to change the scale factor of the sun
planet = 'Sun';
opts.Units = 'mi';
opts.Position = [0, 0, 0];
planet3D(planet, opts);
view([54, 32])


hold on

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
dt = 10000;
tspan = -200*86400:dt:0;
y0 = [r; -VV_minus];
[~, State] = ode113( @(t,y) ode_2body_pb(t,y, mu_Sun), tspan, y0, options );
X_arr = State(:,1); Y_arr = State(:,2); Z_arr = State(:,3);
plot3(X_arr,Y_arr,Z_arr, 'LineWidth', 1.5)

% planet3D('Sun', option);


tspan = 0:dt:200*86400;
y0 = [r; VV_plus];
[~, State] = ode113( @(t,y) ode_2body_pb(t,y, mu_Sun), tspan, y0, options );
X_arr = State(:,1); Y_arr = State(:,2); Z_arr = State(:,3);
plot3(X_arr,Y_arr,Z_arr, 'LineWidth', 1.5)

grid on
axis equal



%% Exercise 1b

mu_Sun = astroConstants(4);
mu_Earth = astroConstants(13);
AU = astroConstants(2);

v_inf_minus = [15.1; 0; 0];
v_inf = norm(v_inf_minus);
Delta = 9200:1000:13200;
r = [1; 0; 0]*AU;

V_Earth = sqrt(mu_Sun/AU);
V_Earth= V_Earth * [0 1 0]';

a = - mu_Earth/v_inf.^2;
turnAng = 2*atan2(-a, Delta);
e = 1./sin(turnAng./2);
r_p = a*(1-e);

delta_V = 2*v_inf*sin(turnAng/2);

VV_minus = V_Earth+ v_inf_minus;


cases = {'leading', 'trailing', 'under'};
chosen_case = cases{1};



u = cross(r, VV_minus);
u = u/norm(u);

    switch chosen_case
        case 'leading'
            rot_dir = -u;
        case 'trailing'
            rot_dir = +u;
        case 'under'
            rot_dir = -V_Earth;      
    end

rot_dir = rot_dir/norm(rot_dir);

v_inf_plus = RodrigueRotationVector(v_inf_minus, turnAng, rot_dir);
VV_plus = V_Earth+ v_inf_plus;

 figure

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
dt = 1000;
tspan = -200*86400:dt:0;
y0 = [r; -VV_minus];
[time, State] = ode113( @(t,y) ode_2body_pb(t,y, mu_Sun), tspan, y0, options );
X_arr = State(:,1); Y_arr = State(:,2); Z_arr = State(:,3);
% plot of the incoming heliptic arc
plot3(X_arr,Y_arr,Z_arr, 'LineWidth', 1.5)
hold on
grid on
axis equal

for i=1:length(Delta)
    tspan = 0:dt:200*86400;
    y0 = [r;VV_plus(:, i)];
    [time, State] = ode113( @(t,y) ode_2body_pb(t,y, mu_Sun), tspan, y0, options );
    X_arr = State(:,1); Y_arr = State(:,2); Z_arr = State(:,3);
    % plot of the outcoming heliptic arc
    plot3(X_arr,Y_arr,Z_arr, 'LineWidth', 1.5)
end



%% Exercise 2
clear all; clc; close all;

mu_Sun = astroConstants(4);
mu_Earth = astroConstants(13);
AU = astroConstants(2);
R_Earth = 6378;
h_atm = 200;

VV_minus = [31.5; 5.2; 0];
VV_plus = [36; 0; 0];

r = [0; -1; 0]*AU;

V_Earth = sqrt(mu_Sun/AU);
V_Earth= V_Earth * [1 0 0]';

vv_inf_minus = VV_minus -V_Earth;
v_inf_minus = norm(vv_inf_minus);

vv_inf_plus =VV_plus -V_Earth;
v_inf_plus = norm(vv_inf_plus);

delta_val = acos(dot(vv_inf_minus, vv_inf_plus)/(v_inf_minus*v_inf_plus));

e = @(r_p, v_inf) 1+ (r_p.*(v_inf).^2)/mu_Earth;

delta = @(r_p, v_inf) 2*asin(1./e(r_p, v_inf));

f_delta = @(r_p)...
            delta_val - (delta(r_p, v_inf_minus)/2 + delta(r_p, v_inf_plus)/2);


r_min = R_Earth + h_atm;
tol = 1e-14;
options = optimoptions('fsolve', 'TolFun', tol, 'TolX', tol);
r_p = fsolve(f_delta, r_min, options);

e_minus = e(r_p, v_inf_minus);
a_minus = r_p/(1-e_minus);

e_plus = e(r_p, v_inf_plus);
a_plus = r_p/(1-e_plus);

delta_v_flyby = vv_inf_plus- vv_inf_minus;

h_GA = r_p - R_Earth;

v_p_minus = sqrt(v_inf_minus^2 + 2*mu_Earth/r_p);
v_p_plus = sqrt(v_inf_plus^2 + 2*mu_Earth/r_p);

delta_V_poweredFB = abs(v_p_plus - v_p_minus);



options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );


% time span
dt = 10; 
tspan = 0:dt:10000;

rr_p = [r_p 0 0]';
vv_p_minus = [0 v_p_minus 0]';
y0 = [rr_p;vv_p_minus];

% integration
[~, State] = ode113( @(t,y) ode_2body_pb(t,y, mu_Earth), -tspan, y0, options );
X = State(:,1); Y = State(:,2); Z = State(:,3);

h = figure;

% plot earth
planet = 'Earth';
opts.Units = 'km';
opts.Position = [0, 0, 0];
planet3D(planet, opts);
view([0,90])

hold on
axis equal
grid on
hold on
plot3(X,Y,Z, 'LineWidth', 2)
X_min_minus = min(X);


rr_p = [r_p 0 0]';
vv_p_plus = [0 v_p_plus 0]';
y0 = [rr_p;vv_p_plus];
[~, State] = ode113( @(t,y) ode_2body_pb(t,y, mu_Earth), tspan, y0, options );
X = State(:,1); Y = State(:,2); Z = State(:,3);
plot3(X,Y,Z, 'LineWidth', 2)
X_min_plus = min(X);



a_minus = abs(a_minus);
a_plus = abs(a_plus);

x0_minus = r_p+a_minus;
x0_plus = r_p+a_plus;
Y_minus = @(x)   tan(pi/2 - delta(r_p, v_inf_minus)/2) * (x-x0_minus);
Y_plus =  @(x)  -tan(pi/2 - delta(r_p, v_inf_plus)/2)  * (x-x0_plus);

X = X_min_minus:x0_minus ;
plot3(X, Y_minus(X), zeros(length(X),1), 'LineWidth', .5, 'Color', 'g')

X =  X_min_plus:x0_plus ;
plot3(X, Y_plus(X), zeros(length(X),1), 'LineWidth', .7, 'Color', 'm')






