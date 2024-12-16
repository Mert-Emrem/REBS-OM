function [a_J2_vect_Margi,a_J2_vect_Ludo] = CHECK_acc_pert_function( kep_el, J2, mu_E, Re)
% function acc_pert_vect computes the accellerationos (perturbations) and
% merge them into a vector 

[r, ~] = par2car(kep_el(1), kep_el(2), kep_el(3), kep_el(4), kep_el(5), kep_el(6), mu_E);

rnorm = norm(r);

z=r(3);


a_J2_vect_Margi = ((3/2* J2 * mu_E * Re^2)/(rnorm^4))*...
       [r(1)/rnorm*(5*(r(3)^2)/(rnorm^2)-1);...
        r(2)/rnorm*(5*(r(3)^2)/(rnorm^2)-1);...
        r(3)/rnorm*(5*(r(3)^2)/(rnorm^2)-3)];

a_J2_vect_Ludo = (3/2) * J2 * (mu_E * Re^2 / rnorm^5) * [(5 * (z/rnorm)^2 - 1) * r(1); (5 * (z/rnorm)^2 - 1) * r(2); (5 * (z/rnorm)^2 - 3) * r(3)];

% % our third body is the sun
% r_CB_3B = uplanet(date+t, 3); % 
% r_SC_3B = 
% a_3B 


end