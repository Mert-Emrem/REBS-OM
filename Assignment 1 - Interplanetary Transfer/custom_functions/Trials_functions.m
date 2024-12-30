%% test functions
% simulator animated plot to see different options

%%
% first attempt
t_dep = 1.506e+04;
ToF1 = 400;
ToF2 = 450;
% Animated_Transfers_Plot([t_dep, t_dep+ToF1, t_dep+ToF1+ToF2])
plotTransfer([t_dep, t_dep+ToF1, t_dep+ToF1+ToF2])
dv = DeltaV_calculator([t_dep, t_dep+ToF1, t_dep+ToF1+ToF2], data, 1)