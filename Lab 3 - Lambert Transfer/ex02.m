clc; clear all;

% Finds total delta v for a Lambert transfer arc

mu_E = astroConstants(13);  

[R_1vec,V_1vec] = kep2car_deg(12500, 0, 0, 0, 0, 120,mu_E);
[R_2vec,V_2vec] = kep2car_deg(9500, 0.3,  0, 0, 0, 250,mu_E);
ToF = 3300;          

muSun = 0.39860e6;      % Sun's gravitational parameter [km^3/s^2];

[a, p, e,ERROR,VI,VF,TPAR,THETA] = lambertMR(R_1vec, R_2vec, ToF, muSun, 0, 0, 0 );


y0 = [R_1vec; V_1vec];

deltaV_1 = vecnorm(VI' - V_1vec);
deltaV_2 = vecnorm(VF' - V_2vec);
deltaV_tot = deltaV_1 + deltaV_2

% Set options
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Set time span
tspan = [0:100:ToF];

% Perform the integration
[T, StateMat] = ode113( @(t,y) ode_2body_pb(t,y,mu_E), tspan, y0, options);

X = StateMat(:,1); 
Y = StateMat(:,2); 
Z = StateMat(:,3);

