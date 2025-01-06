function ds = eq_motion_Gauss_RSW(t, kep_el, acc_pert_function,  mu) %, J2, Re, day)
% ---------------------------------------------------------------------------------------------
% ds = eq_motion_Gauss_RSW(t, kep_el, acc_pert_function,  mu) 
% 
% This function gives as output the system of ODE expressed with the 
% resolution method of Gauss in the
% RSW (radial - transversal - out-of-plane) reference frame, calling the 
% acc_pert_function to calculate the perturbations given by the primary
% planet oblateness and Sun gravitational perturbation as third body
% 
%
% INPUTS
% - ODE variables
%   t                   [1x1]    [s]             time variable of the ODE system 
%   kep_el              [6x1]    [...]           column vector of the 6 keplerian elements variables of the
%                                                ODE system: [a [km]; e [-]; i [rad]; OMEGA [rad]; omega [rad]; theta[rad]]
% - function for acceleration
%   acc_pert_function                            @(t, k_el) acc_pert_function(t, kep_el, J2, mu, R_planet, date_mjd2000, tilt)
%   mu                  [1x1]    [km^3/s^2]      Primary planet gravitational parameter
%
% OUTPUT
%   ds                  [6x1]                    system of 6 ODE in the unknown: a, e, i, OM,
%                                                om, th following the resolution method of Gauss in the
%                                                RSW (radial - transversal - out-of-plane) reference frame
%
% Author: Serlini Mariagiulia, Bernasconi Ludovico, Emrem Mert, Richero
%         Giovanni
% Last update: 23/12/2024
%
% ------------------------------------------------------------------------------------------------

% extract the unknowns from the vector 
a = kep_el(1);
e = kep_el(2); 
i = kep_el(3); 
OM = kep_el(4); 
om = kep_el(5);
th = kep_el(6);

% acceleration vector in cartesian coordinates
a_p = acc_pert_function(t, kep_el);  

% rotational matrix: from cartesian reference frame to RSW
R = iner2RSW(om,th, i, OM);

% acceleration vector expressed in RSW reference frame
a_p_RSW = R*a_p;

a_r = a_p_RSW(1);
a_s = a_p_RSW(2);
a_w = a_p_RSW(3);

% other orbital parameters: semi-latus rectum, magnitude of angular
% momentum, radius magnitude
p = a*(1-e^2);
h = sqrt(p*mu);
r = p/(1+e*cos(th));

% the 6 ODE Gaussian semi-analytical model
da = 2*a^2/h*(e*sin(th)*a_r + p/r * a_s);
de = 1/h * (p*sin(th)*a_r + ((p+r)*cos(th) + r*e)*a_s);
di = r*cos(th+om)/h * a_w;
dOM = r*sin(th+om)/(h*sin(i)) * a_w;
dom = 1/(h*e) * (-p*cos(th)*a_r + (p+r)*sin(th)*a_s) - r*sin(th+om)*cos(i)/(h*sin(i))*a_w;
dth = h/r^2 + 1/(e*h)*(p*cos(th)*a_r - (p+r)*sin(th)*a_s);

% system of 6 ODE
ds = [da; de; di; dOM; dom; dth];

end


