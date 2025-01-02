function [dv] = DeltaV_calculator(x, data, flag_atm)
%
% DeltaV_calculator: Total deltaV of the mission
% 
% Function to compute the total deltaV of the mission, which includes the
% impulsive maneuvers required at departure, flyby, and arrival. 
% It can also consider the atmospheric constraint during the flyby.
% If the flag `flag_atm` is set, the atmospheric constraint is applied to
% the flyby. Otherwise, the pericenter radius of the flyby is unrestricted.
%
% PROTOTYPE:
%  [dv] = DeltaV_calculator(x, data, flag_atm)
% 
% INPUT:
%  x [3,1]        Times array: 
%                    x(1) = departure time [MJD2000]
%                    x(2) = flyby time [MJD2000]
%                    x(3) = arrival time [MJD2000]
%  data           Struct containing mission parameters (e.g., planet radii, 
%                 gravitational parameters, atmospheric conditions)
%  flag_atm       Atmospheric constraint flag:
%                   0 = No atmospheric constraint
%                   1 = Atmospheric constraint applied
% 
% OUTPUT:
%  dv [1]         Total deltaV of the mission [km/s]
%
% AUTHORS: Richero Giovanni, Emrem Mert, Bernarsconi Ludovico, Serlini
% Mariagiulia
%
% -------------------------------------------------------------------------

% Initialize times
t_d = x(1); % Departure time
t_f = x(2); % Flyby time
t_a = x(3); % Arrival time

% Planet IDs
id_d = 1;  % Mercury
id_f = 4;  % Mars
id_a = 40; % Target Asteroid

% Calculate the orbital elements and state vectors of departure, flyby, and arrival bodies
[kep_d, ksun] = uplanet(t_d, id_d); 
[rr_d, vv_d] = kep2car(kep_d(1), kep_d(2), kep_d(3), kep_d(4), kep_d(5), kep_d(6), ksun);

[kep_f, ksun] = uplanet(t_f, id_f); 
[rr_f, vv_f] = kep2car(kep_f(1), kep_f(2), kep_f(3), kep_f(4), kep_f(5), kep_f(6), ksun);

kep_a = ephAsteroids(t_a, id_a); 
[rr_a, vv_a] = kep2car(kep_a(1), kep_a(2), kep_a(3), kep_a(4), kep_a(5), kep_a(6), ksun);

% Solve Lambert's problem for the two transfer legs
[~, ~, e1, err1, vt1_i, vt1_f, ~, ~] = lambertMR(rr_d, rr_f, (t_f - t_d) * 24 * 3600, ksun, 0, 0, 0, 0);
[~, ~, e2, err2, vt2_i, vt2_f, ~, ~] = lambertMR(rr_f, rr_a, (t_a - t_f) * 24 * 3600, ksun, 0, 0, 0, 0);

% Initialize deltaV output
dv = NaN; 

% Ensure there are no errors in the Lambert solutions
if (err1 == 0) && (err2 == 0)
    % Compute the hyperbolic excess velocities (v_infinity) for the flyby
    vinf_m = vt1_f' - vv_f; % Incoming velocity at Mars
    vinf_p = vt2_i' - vv_f; % Outgoing velocity from Mars
    
    % Compute the deltaV required for the powered gravity assist
    if ~flag_atm
        dvp = PowerGravityAssist(vinf_m, vinf_p, data.Mars.Radius, 0, data.Mars.mu, 1);
    else
        dvp = PowerGravityAssist(vinf_m, vinf_p, data.Mars.Radius, data.Mars.h_atm, data.Mars.mu, 1);
        if isnan(dvp)
            disp('!!! Warning: Hyperbolic transfer not feasible due to atmospheric constraint.');
        end
    end
    
    % Compute the total deltaV
    dv = norm(vv_d - vt1_i') + norm(vv_a - vt2_f') + dvp;
end

end
