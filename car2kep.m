function [a, e, i, OM, om, th] = car2kep(rr, vv, mu)
r=norm(rr);
v=norm(vv);
%semiasse maggiore
a=(2/r-v^2/mu)^(-1);
% vettore h
hh=cross(rr,vv);
h=norm(hh);
%vettore eccentricità
ee=cross(vv, hh)./mu-rr./r;
e=norm(ee);
%inclinazione
i = acos(hh(3)/h);
%linea dei nodi
N = cross([0 0 1]', hh)/norm(cross([0 0 1], hh));
if(i==0)
    N=[1 0 0]';
end
%RAAN
OM=acos(N(1))+(N(2)<0)*(2*pi-2*acos(N(1)));
%Anomalia pericentro
om=acos(N'*ee./e)+(ee(3)<0)*(2*pi-2*acos(N'*ee./e));
%Anomalia ffr
th=acos(rr'*ee./(r*e))+(vv'*rr<0)*(2*pi-2*acos(rr'*ee./(r*e)));