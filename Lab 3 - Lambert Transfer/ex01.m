clc;clear all

% Plots the Lambert's arc (prograde, short)

R_1vec = [-21800, 37900, 0];
R_2vec = [27300, 27700, 0];
ToF = 54400;          

mu_E = astroConstants(13);  
muSun = 0.39860e6;      % Sun's gravitational parameter [km^3/s^2];

[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(R_1vec, R_2vec, ToF, muSun, 0, 0, 0 );

R_1vec = R_1vec(:); 
R_2vec = R_2vec(:); 
V_1vec = VI; 
V_2vec = VF;

y0 = [R_1vec; V_1vec'];

% Set options
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Set time span
tspan = [0:100:ToF];

% Perform the integration
[T, StateMat] = ode113( @(t,y) ode_2body_pb(t,y,mu_E), tspan, y0, options);

X = StateMat(:,1); 
Y = StateMat(:,2); 
Z = StateMat(:,3);

figure
plot3(X,Y,Z)
grid on