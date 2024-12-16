function acc_pert_vect = acc_pert_function(t, kep_el, J2, mu_E, Re, date)
% function acc_pert_vect computes the accellerations (perturbations) and
% merge them into a vector 
%
% INPUT:
%  - t                      variable
%  - kep_el     [6x1]       Kelperian elements
%  - J2         [1x1]       Gravitatonal field constant of the Earth 
%  - mu_E       [1x1]       Earth's planetary constant
%  - Re         [1x1]       Earth's radius
%  - date       [1x1]       Julian date of start mission
%
% OUTPUT:
%  - acc_pert_vect  []      Vector containig the disturbance accelerations
%
% Author: Bernasconi Ludovico, Emrem Mert, Richero Giovanni, Serlini
%         Mariagiulia
%
% Date:14/12/2024
%
%--------------------------------------------------------------------------

% Transform orbit from keplerian to cartesian
[r, ~]  = par2car(kep_el(1), kep_el(2), kep_el(3), kep_el(4), kep_el(5), kep_el(6), mu_E);
rnorm   = norm(r);

% Compute J2 perturbation
a_J2_vect   = ((3/2* J2 * mu_E * Re^2)/(rnorm^4))*...
               [r(1)/rnorm*(5*(r(3)^2)/(rnorm^2)-1);...
                r(2)/rnorm*(5*(r(3)^2)/(rnorm^2)-1);...
                r(3)/rnorm*(5*(r(3)^2)/(rnorm^2)-3)];

% % Compute 3 body perturbation with Sun
% r_CB_SC                 = r;                                                                                            %[km]   radius Central Body - Satellite
% [kep_CB_3B , mu_3B]     = uplanet(date+t, 3);                                                                           %[km]   Keplerian Elements Central Body - 3rd Body
% [r_CB_3B,~]             = par2car(kep_CB_3B(1),kep_CB_3B(2),kep_CB_3B(3),kep_CB_3B(4),kep_CB_3B(5),kep_CB_3B(6),mu_3B); %[km]   radius 3rd Body - Central Body
% r_CB_3B                 = - r_CB_3B;                                                                                    %[km]   radius Central Body - Satellite
% r_SC_3B                 = r_CB_3B - r_CB_SC;                                                                            %[km]   radius Satellite - 3rd Body
% 
% c                       = r_SC_3B - r_CB_3B;
% d                       = dot(c,(c - 2*r_SC_3B))/(dot(r_SC_3B,r_SC_3B));
% r_CB_3B_norm            = norm(r_CB_3B);
% d_norm                  = norm(d);
% a_3B                    = mu_3B / r_CB_3B_norm ^3 * ...
%                           (d_norm*(3 + 3*d_norm + d_norm^2)/...
%                           (1 + (1 + d_norm)^(3/2)) ...
%                           * r_SC_3B + c);

% Mux the perturbation in the output vector
acc_pert_vect = a_J2_vect; % + a_3B 
end