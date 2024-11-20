function dy = ode_2body_pb(~, y, mu, J2, Re)

r = y(1:3);
v = y(4:6);
rnorm = norm(r);

% dy = f ( t, y(t) )
dy = [ v;
       (-mu/rnorm.^3)*r];

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