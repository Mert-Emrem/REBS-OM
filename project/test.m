
clear all
close all
clc

addpath ..\LAB
addpath ..\LAB\timeConversion\time\



% 2003 April 1 – 2003 August 1
% dep_date_min = date2mjd2000([2002, 4, 1, 0, 0, 0]);
% dep_date_max = date2mjd2000([2010, 8, 1, 23, 59, 59]);
% 
% dep_span = hms2fracday(23, 0, 0);
% 
% dep_window = dep_date_min: dep_span :dep_date_max;
% 
% constant
% mu_sun = astroConstants(4);
% mu_earth = astroConstants(13);
% mu_mars = astroConstants(14);
% 
% 
% 
% [kep_dep, ~] = uplanet_vec(dep_window, 3);
% 
% [r_dep, v_dep] = kep2car_vec(kep_dep(1,:), kep_dep(2,:), kep_dep(3,:),...
%              kep_dep(4,:), kep_dep(5,:), kep_dep(6,:), mu_sun);
% 
% kep2car(kep_dep(1,1), kep_dep(2,1), kep_dep(3,1),...
%              kep_dep(4,1), kep_dep(5,1), kep_dep(6,1), mu_sun)
% 
% figure
% plot3(r_dep(3, :), r_dep(2, :), r_dep(1, :))
% hold on
% plot3(r_arr_v(1, :), r_arr_v(2, :), r_arr_v(3, :))
% rr = [];
% kep = [];
% 
% for ii = 1:length(dep_window)
%     % Departure Earth
%     [kep_dep, ~] = uplanet(dep_window(ii), 3);
%     kep = [kep, kep_dep'];
% 
%     [r_dep1, v_dep] = kep2car(kep_dep(1), kep_dep(2), kep_dep(3),...
%              kep_dep(4), kep_dep(5), kep_dep(6), mu_sun);
%      rr = [rr, r_dep1]
% 
% end
% 
% a = kep(1,1);
% e = kep(2,1);
% [r_dep, v_dep] = kep2car_vec(a, e, kep(3,:),...
%              kep(4,:), kep(5,:), kep(6,:), mu_sun)
% 
% % modificare matrice mutlidimensionale
% 
% figure
% plot3(r_dep(1,:), r_dep(2, :), r_dep(3, :))
% 
% figure
% plot3(rr(1,:), rr(2, :), rr(3, :))
% 
% %%
% 
% 
% kep2car_vec(1e+6, 0.2, [20 40], [30 50], [20, 30], [10 20], mu_sun)
% kep2car(1e+6, 0.2,  40, 50,  30,  20, mu_sun)
% kep2car(1e+6, 0.2,  20, 30,  20,  20, mu_sun)





%%

% constant
mu_sun = astroConstants(4);
mu_earth = astroConstants(13);
mu_mars = astroConstants(14);


% 2003 April 1 – 2003 August 1
dep_date_min = date2mjd2000([2003, 4, 1, 0, 0, 0]);
dep_date_max = date2mjd2000([2004, 8, 1, 23, 59, 59]);

dep_span = hms2fracday(5, 0, 0);

dep_window = dep_date_min: dep_span :dep_date_max;

% 2003 September 1 – 2004 March 1
arr_time_min = date2mjd2000([2003, 9, 1, 0, 0, 0]);
arr_time_max = date2mjd2000([2004, 3, 1, 23, 59, 59]);

arr_span = hms2fracday(5, 0, 0);

arr_window = arr_time_min: arr_span: arr_time_max;




[kep, ~] = uplanet_vec(dep_window, 3);

a = kep(1,1);
e = kep(2,1);

[r_dep, v_dep] = kep2car_vec(a, e, kep(3,:),...
             kep(4,:), kep(5,:), kep(6,:), mu_sun)

[kep, ~] = uplanet_vec(arr_window, 1);

a = kep(1,1);
e = kep(2,1);

[r_arr, v_arr] = kep2car_vec(a, e, kep(3,:),...
             kep(4,:), kep(5,:), kep(6,:), mu_sun)


figure
plot3(r_dep(1,:), r_dep(2, :), r_dep(3, :))
hold on
plot3(r_arr(1,:), r_arr(2, :), r_arr(3, :))

%%

clear all
close all
clc

% constant
mu_sun = astroConstants(4);
mu_earth = astroConstants(13);
mu_mars = astroConstants(14);


% 2003 April 1 – 2003 August 1
dep_date_min = date2mjd2000([2003, 4, 1, 0, 0, 0]);
dep_date_max = date2mjd2000([2003, 8, 1, 23, 59, 59]);

dep_span = hms2fracday(5, 0, 0);

dep_window = dep_date_min: dep_span :dep_date_max;

% 2003 September 1 – 2004 March 1
arr_time_min = date2mjd2000([2003, 9, 1, 0, 0, 0]);
arr_time_max = date2mjd2000([2004, 3, 1, 23, 59, 59]);

arr_span = hms2fracday(5, 0, 0);

arr_window = arr_time_min: arr_span: arr_time_max;

deltaV = zeros(length(dep_window), length(arr_window));

for ii = 1:length(dep_window)
    % Departure Earth
    [kep_dep, ~] = uplanet(dep_window(ii), 3);

    [r_dep, v_dep] = kep2car(kep_dep(1), kep_dep(2), kep_dep(3),...
             kep_dep(4), kep_dep(5), kep_dep(6), mu_sun);

   
            for jj = 1: length(arr_window)
                % Departure Mars
                [kep_arr, ~] = uplanet(arr_window(jj), 4);

                [r_arr, v_arr] = kep2car(kep_arr(1), kep_arr(2), kep_arr(3),...
                                    kep_arr(4), kep_arr(5), kep_arr(6), mu_sun);
    
                t1 = dep_window(ii);
                t2 = arr_window(jj);

                tof = (t2-t1)*24*3600;
    
                [a_t, p_t, e_t,~, v1_t, v2_t, ~, ~] = lambertMR(r_dep, r_arr, tof, mu_sun, 0, 0, 0, 1);
    
                deltaV(ii,jj) = norm(v1_t - v_dep') + norm(v_arr'-v2_t);
    
            end

        
end

deltaV_opt = min(min(deltaV));
[i, j] = find(deltaV_opt == deltaV);
depart_date = mjd20002date(dep_window(i));
arrival_date = mjd20002date(arr_window(j));

figure
contour(dep_window, arr_window, deltaV', [ 6 7 8 9 10])
hold on
colorbar



