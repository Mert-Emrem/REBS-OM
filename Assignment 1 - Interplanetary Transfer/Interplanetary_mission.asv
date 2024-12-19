clear 
close all
clc

addpath 'functions'\timeConversion\time\
addpath 'functions'\
addpath 'custom_functions'\

%% Interplanetary Explorer Mission

% Gravitational constants of related bodies
mu_sun = astroConstants(4);
mu_merc = astroConstants(11);
mu_mars = astroConstants(14);

R_mars = 3390; % mars radius [km]

%% Departure planet: Mercury

% Synodic period between Mercury-Mars, calculated using Estimation_of_syn_periods.m
T_syn_1 = 100.8882; % [days]

% Departure window in mjd2000
dep_date_min = date2mjd2000([2030, 1, 1, 0, 0, 0]);
dep_date_max = date2mjd2000([2060, 1, 1, 0, 0, 0]);
% dep_date_max = dep_date_min + 10*T_syn_1; % ToF_max

% Time interval for departure window [days]
dep_dt = 50; 

% Create vector of L elements where 
% L = (max. dep. date - min. dep. date)/dt
dep_window = dep_date_min: dep_dt :dep_date_max; % [1 x L]

% Obtain keplerian elements of mercury for each departure date [6 x L]
[kep, ~] = uplanet_vec(dep_window, 1);

% Semimajor axis of mercury from keplerian elements
a_merc = kep(1,1);

% Cartesian states of mercury for each departure date [(3,3) x L]
[r_dep, v_dep] = kep2car_vec(a_merc, kep(2,1), kep(3,1),...
             kep(4,1), kep(5,1), kep(6,:), mu_sun);

%% Flyby planet: Mars

% Flyby window needs not be same as the departure window,
% but it allows for an intuitive visualization
flyby_window = dep_window; % [1 x M]

% Same procedure as Mercury
[kep, ~] = uplanet_vec(flyby_window, 4);

a_mars = kep(1,1);

[r_mars, v_mars] = kep2car_vec(a_mars, kep(2,:), kep(3,:),...
             kep(4,:), kep(5,:), kep(6,:), mu_sun);

%% Arrival asteroid N.40 (Harmonia)

% Synodic period between Mars-Harmonia, calculated using Estimation_of_syn_periods.m
T_syn_2 = 4.1902; % [years]

% Arrival window in mjd2000
arr_time_min = date2mjd2000([2030, 1, 1, 0, 0, 0]);
arr_time_max = date2mjd2000([2060, 1, 1, 0, 0, 0]);

% Time interval for arrival window [days]
arr_dt = 50;

arr_window = arr_time_min: arr_dt: arr_time_max; % [1 x N]

% Obtain keplerian ephemerides of Harmonia
[kep, f, ~] = ephAsteroids_vec(arr_window, 40);

a_harmonia = kep(1);
i = kep(3).*ones(1,length(f));
OM = kep(4).*ones(1, length(f));
om =  kep(5).*ones(1, length(f));

% Cartesian states of Harmonia for each arrival date [(3,3) x L]
[r_harm, v_harm] = kep2car_vec(kep(1), kep(2), i,...
            OM, om, f, mu_sun);

%% First leg of the transfer (Mercury - Mars)

% Preallocation of the matrices to be created
deltaV_Merc_Mars = ones(length(dep_window), length(flyby_window))*NaN;
Vinf_minus = ones(length(dep_window), length(flyby_window),3)*NaN;

for i = 1:length(dep_window)
    for j = 1:length(flyby_window)

        % Compute the time of flight
        ToF = (flyby_window(j) - dep_window(i)) * 86400; % [s]
        
        % Ignore temporally impossible solutions (flyby before departure)
        if ToF > 0
            % Solve Lambert's problem
            [~, ~, ~, err_1, VI, VF, ~, ~] = ...
                lambertMR(r_dep(:, i) , r_mars(:,j) , ...
                ToF, mu_sun, 0, 0, 0);


            % Compute Delta-V
            deltaV_1 = norm(VI' - v_dep(:, i)); % Departure Delta-V
            % deltaV_2 = norm(VF' - v_mars(:,j)); % Arrival Delta-V
            deltaV_2 = 0;

            if err_1 == 0
            deltaV_Merc_Mars(i, j) = deltaV_1+deltaV_2;
            Vinf_minus(i,j,:) = VF' - v_mars(:,j);
            end
        end
    end
end


%% Second leg of the transfer (Mars - Harmonia)

% Preallocation of the matrices to be created
deltaV_Mars_Harm = ones(length(flyby_window), length(arr_window))*NaN;
Vinf_plus = ones(length(flyby_window), length(arr_window),3)*NaN;

for i = 1:length(flyby_window)

    for j = 1:length(arr_window)

        % Compute the time of flight
        ToF = (arr_window(j) - flyby_window(i)) * 86400; % [s]
        
        % Ignore temporally impossible solutions (arrival at Harmonia before Mars flyby)
        if ToF > 0
            % Solve Lambert's problem
            [~, ~, ~, err_2, VI, VF, ~, ~] = ...
                lambertMR(r_mars(:, i) , r_harm(:,j) , ...
                ToF, mu_sun, 0, 0, 0);

            % Compute Delta-V
            deltaV_1 = norm(VI' - v_mars(:,i)); % Departure Delta-V
            deltaV_2 = norm(VF' - v_harm(:,j)); % Arrival Delta-V
            deltaV_1 = 0;

            if err_2 == 0
               deltaV_Mars_Harm(i, j) = deltaV_1+deltaV_2;
               Vinf_plus(i,j,:) = VI' - v_mars(:,i);
            end
        end
    end

end

%% Post-processing of the first and seconds legs, Gravity Assist

% Repeat DeltaV results of the first leg over every arrival date
Merc_Mars_3d = repmat(deltaV_Merc_Mars, [1, 1, length(arr_window)]);

% Repeat DeltaV results of the second leg over every departure date
Mars_Harm_3d = repmat(deltaV_Mars_Harm, [1, 1, length(dep_window)]);

% Reorient the DeltaV matrix of first leg such that the two are combinable
Merc_Mars_3d = permute(Merc_Mars_3d, [2, 3, 1]);

% Generate 2D and 3D porkchop plots of both legs of the transfer
porkchopPlotter(deltaV_Merc_Mars, flyby_window, dep_window)
porkchopPlotter(deltaV_Mars_Harm, arr_window, dep_window)

% Briefer names for ease of use
M1 = Mars_Harm_3d;
M2 = Merc_Mars_3d;

%% Optimization

Delta_GA = ones(length(flyby_window), length(arr_window), length(dep_window))*NaN;
DeltaVtot = ones(length(flyby_window), length(arr_window), length(dep_window))*NaN;

for i=1:length(flyby_window)
    for j=1:length(arr_window)
        for k=1:length(dep_window)

             if dep_window(k)<flyby_window(i) && flyby_window(i)<arr_window(j)

             vinf_m = squeeze(Vinf_minus(k,i,:));
             vinf_p = squeeze(Vinf_plus(i,j,:));
             dvp  = PowerGravityAssist(vinf_m, vinf_p...
                    ,R_mars, 100, mu_mars);

             Delta_GA(i, j, k) = dvp;

                  if not(isnan(dvp)) &&...
                     not(isnan(M1(i, j, k))) && ...
                     not(isnan(M2(i, j, k)))
    
                     DeltaVtot(i,j,k) = dvp + M1(i, j, k) + M2(i, j, k);
    
                  end

             end

        end

    end

end


%%

DeltaV_3dofs_Plotter(DeltaVtot, 200, 180)
[Opt, idx] = min(DeltaVtot(:));
[row, col, depth] = ind2sub(size(DeltaVtot), idx);
t_flyby = flyby_window(row);
t_arr = arr_window(col);
t_dep = dep_window(depth);


%% plot

plotTransfer([t_dep, t_flyby, t_arr],r_dep, r_mars, r_harm)

%
figure;
% 
planet = 'Sun';
opts.Units = 'km';
opts.Position = [0, 0, 0];

planet3D(planet, opts);
view([54, 32])
hold on

x = r_dep(1,1);
y = r_dep(2,1);
z = r_dep(3,1);

% line=animatedline(x, y, z,'color', '#ffff00', 'LineWidth',2);
line = animatedline;

for t = 1:size(r_dep,2)

        x = r_dep(1, t);
        y = r_dep(2, t);
        z = r_dep(3, t);

        addpoints(line, x, y, z);

        drawnow;

        pause(0.01);

end


%%

N_runs = 5; % Number of runs
% To check algorithm convergence set N_runs > 1.

data.R_mars = R_mars;
data.mu_mars = mu_mars;

% Save results for each run:
x_runs = [];
dv_runs = [];

disp(['fmincon search with ',num2str(N_runs),' runs running..']);
tic

x0 = [t_dep, t_flyby-t_dep, t_arr-t_flyby];


data.h = 100;

for i = 1:N_runs

options=optimoptions("fmincon");
[x, dv] = fmincon(@(x) dvFun(x, data), x0, [], [], [], [],...
            [dep_window(1), 10, 10],...
            [dep_window(end), 1000,10000],...
            @(x) nonlcon(x, data), options);
x = cumsum(x);

dv_runs = [dv; dv];
x_runs = [x_runs; x];

end

depdate = mjd20002date(x(1));
fprintf( ['Optimized departure date from Mercury: ', repmat('%d ', 1, numel(depdate)), '\n'], depdate);

flybydate = mjd20002date(x(2));
fprintf( ['Optimized flyby date via Mars: ', repmat('%d ', 1, numel(depdate)), '\n'], flybydate);

arrdate = mjd20002date(x(3));
fprintf( ['Optimized arrival date to Harmonia: ', repmat('%d ', 1, numel(arrdate)), '\n'], arrdate);


toc

% Select minimum deltaV solution:

[~,index] = min(dv_runs);
dv = dv_runs(index)
x = x_runs(index,:);

plotTransfer([x(1), x(2), x(3)],r_dep, r_mars, r_harm)

%%

N_runs = 5; % Number of runs
% To check algorithm convergence set N_runs > 1.

data.R_mars = R_mars;
data.mu_mars = mu_mars;

% Save results for each run:
x_runs = [];
dv_runs = [];

% disp(['ga search with ',num2str(N_runs),' runs running..']);
% tic
% 
% data.h = 100;
% A = [1 1 1];
% b = arr_window(end);
% 
% for i = 1:N_runs
% 
% options=optimoptions("ga");
% [x_ga, dv] = ga(@(x) dvFun(x, data), 3,...
%                 A, b, [], [],...
%                 [dep_window(1), 10, 10],...
%                 [dep_window(end), 1000,10000],...
%                 @(x) nonlcon(x, data), options);
%                 x = cumsum(x);
%             % [x(1)-500, x(2)-500-x(1), x(3)-1000-x(2)],...
%             % [x(1)+500, x(2)+500-x(1), x(3)+1000-x(2)]);
%             % % @(x) nonlcon(x, data), options);
% 
% x_ga = cumsum(x_ga);
% 
% dv_runs = [dv; dv];
% x_runs = [x_runs; x_ga];
% 
% end
% 
% toc


% Select minimum deltaV solution:

[~,index] = min(dv_runs);
dv = dv_runs(index)
x = x_runs(index,:);


% Vinf_minus_val = vecnorm(Vinf_minus,2,1);
% Vinf_plus_val = vecnorm(Vinf_plus,2,1);
% 
% delta_val = acos(dot(Vinf_minus, Vinf_plus)./(Vinf_minus_val.*Vinf_plus_val));
% 
% e = @(r_p, v_inf) 1+ (r_p.*(v_inf).^2)/mu_mars;
% 
% delta = @(r_p, v_inf) 2*asin(1./e(r_p, v_inf));
% 
% f_delta = @(r_p)...
%             delta_val - (delta(r_p, Vinf_minus_val)/2 + delta(r_p, Vinf_plus_val)/2);
% 
% h_atm = 100000;
% 
% r_min = (R_mars + h_atm).*ones(1,length(flyby_window));
% 
% if (f_delta(r_min)<=0)&&(f(1e10)>0)
% 
% tol = 1e-14;
% options = optimoptions('fsolve', 'TolFun', tol, 'TolX', tol);
% 
% r_p = fsolve(f_delta, r_min, options);
% 
% 
% e_minus = e(r_p, Vinf_minus_val);
% a_minus = r_p/(1-e_minus);
% 
% e_plus = e(r_p, Vinf_plus_val);
% a_plus = r_p/(1-e_plus);
% 
% delta_v_flyby = Vinf_minus_val- Vinf_plus_val;
% 
% h_GA = r_p - R_mars;
% 
% v_p_minus = sqrt(Vinf_minus_val.^2 + 2*mu_mars./r_p);
% v_p_plus = sqrt(Vinf_plus_val.^2 + 2*mu_mars./r_p);
% 
% delta_V_poweredFB = abs(v_p_plus - v_p_minus);
% 
% else
%     delta_V_poweredFB = NaN;
%     r_p = NaN;
% 
% end

