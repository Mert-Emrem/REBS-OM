function acc_pert_vect = acc_pert_function_cartesian(t, cart_coord, J2, mu, R_planet, date_mjd2000, tilt)

% perturbation vector in cartesian coordinates: a_J2 
r = cart_coord(1:3);

rnorm = norm(r);


a_J2_vect = ((3/2* J2 * mu * R_planet^2)/(rnorm^4))*...
       [r(1)/rnorm*(5*(r(3)^2)/(rnorm^2)-1);...
        r(2)/rnorm*(5*(r(3)^2)/(rnorm^2)-1);...
        r(3)/rnorm*(5*(r(3)^2)/(rnorm^2)-3)];

% the third body is the Sun
% central planet position wrt the Sun in the ecliptic RF
TIME = date_mjd2000 + t/(60*60*24);
[kep,mu_3B] = uplanet(TIME, 1);
[r_3B_CB, ~]  = par2car(kep(1), kep(2), kep(3), kep(4), kep(5), kep(6), mu_3B); % cartesian centered in the sun
%r_CB_3B = -r_CB_3B;

% position of the Sun wrt the central planet
rot_frame = [1, 0, 0; 0 cos(tilt) sin(tilt); 0, -sin(tilt) cos(tilt)];
r_CB_3B = rot_frame*(-r_3B_CB);

% S/C position wrt the central planet
r_CB_SC = r; % cartesian coordinates centered in the earth

% S/C position wrt the Sun
r_SC_3B = r_CB_3B - r_CB_SC;

% perturbation
c = r_SC_3B - r_CB_3B;
d_vect = dot(c, c-2*r_SC_3B)/dot(r_SC_3B, r_SC_3B);
d = norm(d_vect);
diff = 1/norm(r_CB_3B)^3 *(d*(3+3*d+d^2)/(1+(1+d)^1.5)*r_SC_3B + c);

a_3B = mu_3B *diff;

acc_pert_vect = a_J2_vect + a_3B;
end