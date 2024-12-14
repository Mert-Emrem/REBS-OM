function dy = ode_2body_pb(~, y, mu, J2, Re)
%
% this function creates the system of ode to solve the 2 body problem
%
% the fucntion works for both non perturbed and perturbed orbits
% 
% NON PERTURBED
% 
%     INPUT
%     t[1] Time (can be omitted, as the system is autonomous) [T]
%     y[6x1] State of the body ( rx, ry, rz, vx, vy, vz ) [ L, L/T ]
%     mu[1] Gravitational parameter of the primary [L^3/T^2]
%     
%     OUTPUT:
%     dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]
% 
% PERTURBED
%
%     ADDITIONAL INPUT
%     J2 --> constant for pertubations
%     Re --> Radius of the attractor
%     
%     OUTPUT:
%     dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]


r = y(1:3);
v = y(4:6);
rnorm = norm(r);

% non perturbed
dy = [ v;
       (-mu/rnorm.^3)*r];

% perturbed
if nargin>3
a_J2 = ((3/2* J2 * mu * Re^2)/(rnorm^4))*...
       [r(1)/rnorm*(5*(r(3)^2)/(rnorm^2)-1);...
        r(2)/rnorm*(5*(r(3)^2)/(rnorm^2)-1);...
        r(3)/rnorm*(5*(r(3)^2)/(rnorm^2)-3)];

dy = [ v(1);
       v(2);
       v(3);
       (-mu/rnorm^3)*r(1) + a_J2(1);
       (-mu/rnorm^3)*r(2) + a_J2(2);
       (-mu/rnorm^3)*r(3) + a_J2(3)];

end