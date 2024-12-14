%
%% Kepler's equation

a = 7000*1e+3;
mu_earth = astroConstants(13);
k = 2;
N = 100;
e = [0, 0.2, 0.4, 0.6, 0.8, 0.95];
% e = 0:0.01:0.98;
e = e';
t0 = 0;
E0 = 0;

n = sqrt(mu_earth/a^3);
T = 2*pi/n;

t = linspace(t0, 2*T, N);

E = KeplerSolver(t, e, a, mu_earth, t0, E0, 1e-4);

E = E/(3pi)*180;


% figure
% plot(t, E)
% grid on;

t_n = t./T;

figure
plot(t_n, E(1, :), 'LineStyle','-')
hold on
plot(t_n, E(2, :), 'LineStyle', '-')
plot(t_n, E(3, :), 'LineStyle', '--')
plot(t_n, E(4, :), 'LineStyle', '-')
plot(t_n, E(5, :))
plot(t_n, E(6, :), 'LineStyle', '-.')
grid on
xlabel('t [T]')
ylabel('E [deg]')
legend('e=0', 'e=0.2', 'e=0.4', 'e=0.6', 'e=0.8', 'e=0.95')

figure
surf(t, e, E);