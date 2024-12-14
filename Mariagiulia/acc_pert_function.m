function acc_pert_vect = acc_pert_function(t, kep_el, J2, mu_E, Re, date_mjd2000)

% perturbation vector in cartesian coordinates: a_J2 
[r, ~] = par2car(kep_el(1), kep_el(2), kep_el(3), kep_el(4), kep_el(5), kep_el(6), mu_E);

rnorm = norm(r);


a_J2_vect = ((3/2* J2 * mu_E * Re^2)/(rnorm^4))*...
       [r(1)/rnorm*(5*(r(3)^2)/(rnorm^2)-1);...
        r(2)/rnorm*(5*(r(3)^2)/(rnorm^2)-1);...
        r(3)/rnorm*(5*(r(3)^2)/(rnorm^2)-3)];
% our third body is the sun

% [kep,mu_3B] = uplanet(date_mjd2000 + t, 3);
% [r_CB_3B, ~]  = par2car(kep(1), kep(2), kep(3), kep(4), kep(5), kep(6), mu_3B); % cartesian centered in the sun
% r_CB_3B = -r_CB_3B;
% r_CB_SC = r; % cartesian coordinates centered in the earth
% r_SC_3B = r_CB_3B - r_CB_SC;
% 
% c = r_SC_3B - r_CB_3B;
% d_vect = dot(c, c-2*r_SC_3B)/dot(r_SC_3B, r_SC_3B);
% d = norm(d_vect);
% diff = 1/norm(r_CB_3B)^3 *(d*(3+3*d+d^2)/(1+(1+d)^1.5)*r_SC_3B + c);
% 
% a_3B = mu_3B *diff;

acc_pert_vect = a_J2_vect; % + a_3B;
end