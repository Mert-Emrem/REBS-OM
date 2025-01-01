function a = a_groundTrack(k, m, omega_Me)

% k = satellite revolution
% m = earth revolution

mu = astroConstants(11);
%omega_Me = 15.04/180 * pi/60/60;

n = omega_Me*k/m;
a = (mu/n^2)^(1/3);


end


