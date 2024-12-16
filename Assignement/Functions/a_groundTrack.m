function a = a_groundTrack(k, m, mu_planet, OMEGA_planet)

n = OMEGA_planet*k/m;
a = (mu_planet/n^2)^(1/3);

end


