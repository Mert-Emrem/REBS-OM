function [c, ceq] = nonlcon(x, data)
%
% nonlcon
%
% PROTOTYPE:
% [c, ceq] = nonlcon(x, data)
%
% INPUT:
%  x [3,1]        Times array:         x(1) = departure time [MJD2000]
%                                      x(2) = flyby time [MJD2000]
%                                      x(3) = arrival time [MJD2000]
%
% OUTPUT:
%  c              Array of nonlinear inequalities (c < 0)
%  ceq            Array of nonlinear equalities (ceq = 0)
%
% AUTHORS: Bernasconi Ludovico, Emrem Mert, Richero Giovanni, Serlini
%          Mariagiulia
%
% -------------------------------------------------------------------------


% Initialize times:
t_d = x(1); t_f = x(2); t_a = x(3); 

% Planets id numbers:
id_d = 1;                              % Mercury
id_f = 4;                              % Mars
id_a = 40;                             % Asteroid

[kep_d, ksun] = uplanet(t_d, id_d);
[rr_d, ~] = kep2car(kep_d(1), kep_d(2), kep_d(3), kep_d(4), kep_d(5), kep_d(6), ksun);
[kep_f, ksun] = uplanet(t_f, id_f);
[rr_f, vv_f] = kep2car(kep_f(1), kep_f(2), kep_f(3), kep_f(4), kep_f(5), kep_f(6), ksun);
kep_a = ephAsteroids(t_a, id_a);
[rr_a, ~] = kep2car(kep_a(1), kep_a(2), kep_a(3), kep_a(4), kep_a(5), kep_a(6), ksun);

[~, ~, e1, err1, ~, vt1_f, ~, ~] = lambertMR(rr_d, rr_f, (t_f - t_d) * 24 * 3600, ksun, 0, 0, 0, 0);
[~, ~, e2, err2, vt2_i, ~, ~, ~] = lambertMR(rr_f, rr_a, (t_a - t_f) * 24 * 3600, ksun, 0, 0, 0, 0);

rp = 0; % Initialize output

% Check: there must be no errors in lambertMR
if (err1 == 0) && (err2 == 0)
    vinf_m = vt1_f' - vv_f;
    vinf_p = vt2_i' - vv_f;
    [~, ~, rp] = PowerGravityAssist(vinf_m, vinf_p, ...
                    data.Mars.Radius, data.Mars.h_atm, data.Mars.mu, 1);
end

if isnan(rp) || ~isfinite(rp) || ~isreal(rp)
    rp = 0;
end

% Nonlinear equality constraints
ceq = 0;

% Nonlinear inequality constraints
c1 = -e1;                  % Ensure e1 >= 0
c2 = e1 - 0.9999;          % Ensure e1 < 1
c3 = -e2;                  % Ensure e2 >= 0
c4 = e2 - 0.9999;          % Ensure e2 < 1
c5 = data.Mars.Radius + data.Mars.h_atm - rp + 1e-2; % Additional rp constraint

% Combine all inequality constraints
c = [c1; c2; c3; c4; c5];
