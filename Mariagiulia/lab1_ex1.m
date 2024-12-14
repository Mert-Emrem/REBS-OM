
mu_E = astroConstants(13);

%% exercise 1

% 1
r0_1 = [26578.137; 0; 0];
v0_1 = [0; 2.221; 3.173];
y0_1 = [r0_1, v0_1];

OrbitAnalysis(r0_1, v0_1, y0_1, mu_E)

%% 2
r0_2 = [6495, -970, -3622 ];
v0_2 = [4.752,2.130,7.950];
y0_2 = [r0_2, v0_2];

OrbitAnalysis(r0_2, v0_2, y0_2, mu_E)
