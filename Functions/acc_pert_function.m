function acc_pert_vect = acc_pert_function(t, kep_el, J2, mu_E, Re)
% function acc_pert_function computes the accellerationos (perturbations) and
% merge them into a vector 
%
% INPUT:
% - t: time
% - kep_el: keplerian elements
% - 
[r, ~] = par2car(kep_el(1), kep_el(2), kep_el(3), kep_el(4), kep_el(5), kep_el(6), mu_E);

rnorm = norm(r);


a_J2_vect = ((3/2* J2 * mu_E * Re^2)/(rnorm^4))*...
       [r(1)/rnorm*(5*(r(3)^2)/(rnorm^2)-1);...
        r(2)/rnorm*(5*(r(3)^2)/(rnorm^2)-1);...
        r(3)/rnorm*(5*(r(3)^2)/(rnorm^2)-3)];

% % our third body is the sun
% r_CB_3B = uplanet(date+t, 3); % 
% r_SC_3B = 
% a_3B 


acc_pert_vect = a_J2_vect; % + a_3B;
end