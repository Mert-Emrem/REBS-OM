function acc_pert_vect = acc_pert_function_cartesian_test_value(t, cart_coord, J2, mu, R_planet, date_mjd2000, tilt, n_planet)
% ---------------------------------------------------------------------------------------------
% acc_pert_vect = acc_pert_function_cartesian(t, cart_coord, J2, mu, R_planet, date_mjd2000, tilt)
% 
% This function gives as output the vector of the perturbating acceleration
% due to primary planet oblateness (J2 effect) and Sun perturbation as third body 
% The function is called inside the function eq_motion_cartesian that
% creates the ode for the semi-analitical resolution method of Gauss in the
% cartesian reference frame
%
% INPUTS
% - ODE variables
%   t               [1x1]    time variable of the ODE system 
%   cart_coord      [6x1]    column vector of radius and velocity vectors variables of the
%                            ODE system: [rx [km]; ry [km]; rz [km]; vx [km/s]; vy [km/s]; vz[km/s]]
% - parameters
%   J2              [1x1]    Gravitational Harmonic coefficent of the primary planet
%   mu              [1x1]    Primary planet gravitational parameter
%   R_planet        [1x1]    primary planet radius
%   date_mjd200     [1x1]    starting date of the orbit evaluated in mjd2000
%   tilt            [1x1]    tilt angle of the primary planet with respect to the ecliptic plane
%  n_planet         [1x1]    number of the primary planet compatible with uplanet function
%
% OUTPUT
%   acc_pert_vect   [3x1]    vector of the perturbing acceleration, sum of the contribution given 
%                            by J2 and the 3BP perturbation of the Sun in
%                            the cartesian reference frame
%
% Author: Serlini Mariagiulia
% Last update: 27/12/2024
%
% ------------------------------------------------------------------------------------------------


% perturbation vector in cartesian coordinates: a_J2 
r = cart_coord(:, 1:3);

rnorm = vecnorm(r);

a_J2_vect = ((3/2* J2 * mu * R_planet^2)./(rnorm.^4))*...
       [r(1, :)./rnorm.*(5*(r(3, :).^2)./(rnorm.^2)-1);...
        r(2, :)./rnorm.*(5*(r(3, :).^2)./(rnorm.^2)-1);...
        r(3, :)./rnorm.*(5*(r(3, :).^2)./(rnorm.^2)-3)];

% the third body is the Sun
% central planet position wrt the Sun in the ecliptic RF
TIME = date_mjd2000+ t./86400;
for i = 1:size(TIME)
    [kep_v, mu_3B] = uplanet(TIME(i), n_planet);
    kep(1,:) = kep_v;
end
for i = 1:size(kep, 1)
    [r_3B_CB_v, ~]  = par2car(kep(i, 1), kep(i, 2), kep(i, 3), kep(i, 4), kep(i, 5), kep(i, 6), mu_3B); % cartesian centered in the sun
    r_3B_CB(i,:) = r_3B_CB_v;
end
% position of the Sun wrt the central planet
rot_frame = [1, 0, 0; 0 cos(tilt) sin(tilt); 0, -sin(tilt) cos(tilt)];
for i = 1:size(r_3B_CB, 1)
    r_CB_3B_v = rot_frame*(-r_3B_CB(i, :))';
    r_CB_3B(i, :) = r_CB_3B_v;
end


% S/C position wrt the central planet
r_CB_SC = r; % cartesian coordinates centered in mercury

% S/C position wrt the Sun
r_SC_3B = r_CB_3B - r_CB_SC;

% perturbation
c = r_SC_3B - r_CB_3B;
d_vect = dot(c, c-2.*r_SC_3B)/dot(r_SC_3B, r_SC_3B);
d = vecnorm(d_vect);
diff = 1./vecnorm(r_CB_3B)^3 .*(d.*(3+3.*d+d.^2)/(1+(1+d).^1.5).*r_SC_3B + c);

a_3B = mu_3B .*diff;

acc_pert_vect = a_J2_vect + a_3B;

end

% a_J2_vect = ((3/2* J2 * mu * R_planet^2)/(rnorm^4))*...
%        [r(1)/rnorm*(5*(r(3)^2)/(rnorm^2)-1);...
%         r(2)/rnorm*(5*(r(3)^2)/(rnorm^2)-1);...
%         r(3)/rnorm*(5*(r(3)^2)/(rnorm^2)-3)];
% 
% % the third body is the Sun
% % central planet position wrt the Sun in the ecliptic RF
% TIME = date_mjd2000 + t/(60*60*24);
% [kep,mu_3B] = uplanet(TIME, n_planet);
% [r_3B_CB, ~]  = par2car(kep(1), kep(2), kep(3), kep(4), kep(5), kep(6), mu_3B); % cartesian centered in the sun
% 
% % position of the Sun wrt the central planet
% rot_frame = [1, 0, 0; 0 cos(tilt) sin(tilt); 0, -sin(tilt) cos(tilt)];
% r_CB_3B = rot_frame*(-r_3B_CB);
% 
% % S/C position wrt the central planet
% r_CB_SC = r; % cartesian coordinates centered in the earth
% 
% % S/C position wrt the Sun
% r_SC_3B = r_CB_3B - r_CB_SC;
% 
% % perturbation
% c = r_SC_3B - r_CB_3B;
% d_vect = dot(c, c-2*r_SC_3B)/dot(r_SC_3B, r_SC_3B);
% d = norm(d_vect);
% diff = 1/norm(r_CB_3B)^3 *(d*(3+3*d+d^2)/(1+(1+d)^1.5)*r_SC_3B + c);
% 
% a_3B = mu_3B *diff;
% 
% acc_pert_vect = a_J2_vect + a_3B;
% end