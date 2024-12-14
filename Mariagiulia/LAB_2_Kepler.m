%% LAB 2
% Kepler equation to calculate E starting from 

a = 7000e+3; 
mu_E = astroConstants(13);
E0 = 0; %[rad]
k = 2;
N = 100;
e = [0; 0.2; 0.4; 0.6; 0.8; 0.95];

n = sqrt(mu_E/a^3);
T = 2*pi/n;

t0 = 0;
t = linspace(t0, 2*T, N); 

figure
grid on 
hold on
title ('E(t,e)');
xlabel('Time [s]');
ylabel('E [deg]');
E_mat = zeros(length(e), length(t));
for i = 1:length(e)
    E = KeplerFunction(t, e(i), a, mu_E, t0, E0, 1e-4);
    E = E/pi*180; %[from rad to deg]
    plot(t, E, '--')
    hold on
    E_mat(i, :) = E;
end
legend ('e=0', 'e=0.2', 'e=0.4', 'e=0.6', 'e=0.8', 'e=0.95');

figure
surf(t, e, E_mat)

%%
E = KeplerFunction(t, e, a, mu_E, t0, E0, 1e-4);


