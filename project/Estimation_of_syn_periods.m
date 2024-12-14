%% synodic period

mu_sun = astroConstants(4);
a_merc = 5.7909e+07;
a_mars = 2.2794e+08;
a_harmonia = 3.3922e+08;

% synodic period of mercury-asteroid
T1 = 2*pi*sqrt(a_merc^3/mu_sun);
T2 = 2*pi*sqrt(a_mars^3/mu_sun);
T_syn_1 = T1*T2/abs(T1-T2);
T_syn_1 = T_syn_1/3600/24;

% tof1 =

% synodic period of asteroid-mars
T2 = 2*pi*sqrt(a_mars^3/mu_sun);
T3 = 2*pi*sqrt(a_harmonia^3/mu_sun);

T_syn_2 = T3*T2/abs(T3-T2);
T_syn_2 = T_syn_2/3600/24/365;

% tof2 =

%% Hohmann


v_mars = sqrt(mu_sun/a_mars);
v_harmonia = sqrt(mu_sun/a_harmonia);

a_h_1 = (a_mars + a_harmonia)/2;

% delta_v1 = abs(v_mars - sqrt(2/a_mars - 1/a_h));
% delta_v2 = abs(v_harmonia - sqrt(2/a_harmonia - 1/a_h));
% 
% deltaVtot = delta_v1 + delta_v2;

T = pi*sqrt(a_h_1^3/mu_sun);

T_max = 1.2*T_syn_2;

%%


v_mars = sqrt(mu_sun/a_mars);
v_merc = sqrt(mu_sun/a_merc);

a_h_2 = (a_mars + a_merc)/2;

% delta_v1 = abs(v_mars - sqrt(2/a_mars - 1/a_h));
% delta_v2 = abs(v_mercury - sqrt(2/a_mercury - 1/a_h));
% deltaVtot = delta_v1 + delta_v2;

T = pi*sqrt(a_h_2^3/mu_sun);