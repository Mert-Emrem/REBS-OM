
mu_E = astroConstants(13);
Re = astroConstants(23);
J2 = astroConstants(9);

%% exercise 2

% 1
r0_1 = [26578.137; 0; 0];
v0_1 = [0; 2.221; 3.173];
y0_1 = [r0_1, v0_1];

% 5 periodi
OrbitAnalysis(r0_1, v0_1, y0_1, mu_E, J2, Re)

% 1 anno 
% 365×24×60×60=31.536.000
% tspan = 0:1000:31536000;
% tspan = linspace( 0, 31536000, 731*1000 );
tspan = linspace( 0, 31536000, 10000 );
OrbitAnalysis(r0_1, v0_1, y0_1, mu_E, J2, Re, tspan)


%% 2
r0_2 = [6495, -970, -3622 ];
v0_2 = [4.752,2.130,7.950];
y0_2 = [r0_2, v0_2];

% 5 periodi
OrbitAnalysis(r0_2, v0_2, y0_2, mu_E, J2, Re)
% 
% % 1 anno 
% % 365×24×60×60=31.536.000
% % tspan = 0:1000:31536000;
% 
% % tspan = linspace( 0, 31536000, 731*1000 );
tspan = linspace( 0, 31536000, 10000 );
OrbitAnalysis(r0_2, v0_2, y0_2, mu_E, J2, Re, tspan)



