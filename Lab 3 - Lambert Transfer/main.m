

% constant
mu_sun = astroConstants(4);
mu_mars = astroConstants(14);

% mean radius sun-mercury
R1 = 57.9e+6;
T1 = 2*pi*sqrt(R1^3/mu_sun);
T1 = T1/3600/24;

% mean radius sun-mars
R2 = 227.9e+6;
T2 = 2*pi*sqrt(R2^3/mu_sun);
T2 = T2/3600/24;

% synodic period of mercury-asteroid
T1 = 2*pi*sqrt(R1^3/mu_sun);
T2 = 2*pi*sqrt(R2^3/mu_sun);
T_syn_1 = abs(1/T1 - 1/T2)

% tof1 =

% synodic period of asteroid-mars
T2 = 2*pi*sqrt(R2^3/mu_sun);
T3 = 2*pi/mean(M);
T_syn_2 = abs(1/T3 - 1/T2)