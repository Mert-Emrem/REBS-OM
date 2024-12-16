function a = a_groundTrack(k, m)

% k = satellite revolution
% m = earth revolution

mu = astroConstants(13);
omegaE = 15.04/180 * pi/60/60;

n = omegaE*k/m;
a = (mu/n^2)^(1/3);

end


