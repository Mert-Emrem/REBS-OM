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

data.Mars.Radius = 3390; % mars radius [km]
data.Mars.mu = mu_mars;
data.Mars.h_atm = 0;

%% Departure planet: Mercury

% Synodic period between Mercury-Mars, calculated using Estimation_of_syn_periods.m
T_syn_1 = 100.8882; % [days]

% Departure window in mjd2000
dep_date_min = date2mjd2000([2030, 1, 1, 0, 0, 0]);
dep_date_max = date2mjd2000([2044, 4, 1, 0, 0, 0]);
% dep_date_max = dep_date_min + 10*T_syn_1; % ToF_max

% Time interval for departure window [days]
dep_dt = 15; 

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
arr_time_max = date2mjd2000([2044, 4, 1, 0, 0, 0]);

% Time interval for arrival window [days]
arr_dt = 15;

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

            if err_1 ~= 1 && err_1 ~= 3 && err_1 ~= 4
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

            if err_2 ~= 1 && err_2 ~= 3 && err_2 ~= 4
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

% Briefer names for ease of use
M1 = Mars_Harm_3d;
M2 = Merc_Mars_3d;

%%

porkchopPlotter1(deltaV_Merc_Mars, flyby_window, dep_window)
porkchopPlotter2(deltaV_Mars_Harm, arr_window, dep_window)

%%
% Generate 2D and 3D porkchop plots of both legs of the transfer

DeltaV_3dofs_Plotter(deltaV_Merc_Mars, deltaV_Mars_Harm, dep_window, flyby_window, arr_window)

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
            ,data.Mars.Radius, data.Mars.h_atm, data.Mars.mu);

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

% DeltaV_3dofs_Plotter(DeltaVtot, 1000, 180)
[Opt, idx] = min(DeltaVtot(:));
[row, col, depth] = ind2sub(size(DeltaVtot), idx);
t_flyby = flyby_window(row);
t_arr = arr_window(col);
t_dep = dep_window(depth);


%% plot

plotTransfer([t_dep, t_flyby, t_arr])


%% animated plot

Animated_Transfers_Plot([t_dep, t_flyby, t_arr])


%%

NNtrials = 5; % Number of trials
% To check algorithm convergence set NNtrials > 1.


% Save results for each run:
x_trials = zeros([NNtrials, 3]);
dv_trials = zeros([NNtrials, 1]);

disp(['fmincon search with ',num2str(NNtrials),' trials running..']);
tic

x0 = [t_dep, t_flyby, t_arr];


for i = 1:NNtrials

options=optimoptions("fmincon");
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

plotTransfer([x(1), x(2), x(3)])

%%

NNtrials = 5; % Number of trials
% To check algorithm convergence set NNtrials > 1.


% Save results for each run:
x_trials = [];
dv_trials = [];

% disp(['ga search with ',num2str(NNtrials),' trials running..']);
% tic
% 
% data.h = 100;
% A = [1 1 1];
% b = arr_window(end);
% 
% for i = 1:NNtrials
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
% dv_trials = [dv; dv];
% x_trials = [x_trials; x_ga];
% 
% end
% 
% toc


% Select minimum deltaV solution:

[~,index] = min(dv_trials);
dv = dv_trials(index)
x = x_trials(index,:);







