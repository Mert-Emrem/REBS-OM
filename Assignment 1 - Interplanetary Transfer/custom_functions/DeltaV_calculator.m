function [dv] = DeltaV_calculator(x, data, flag_atm)
%
% DeltaV_calculator: Total deltaV of the mission
% 
% Function to compute the total deltaV of the mission given the times of
% departure, flyby and arrival. In this case, there is no
% constraint on the minimum value of radius of pericenter of the flyby
% manouvre (rp_min = 0).
%
% PROTOTYPE:
%  [dv, depdate] = dvFun(x, data)
% 
% INPUT:
%  x [3,1]        Times array:         x(1) = departure time [MJD2000]
%                                      x(2) = flyby time [MJD2000]
%                                      x(3) = arrival time [MJD2000]
%  data           Struct containing mission parameters for the flyby
%                 trajectory
%
%  flag_atm           if flag_atm == 1 constraint on 
% 
% OUTPUT:
%  dv [1]         Total deltaV of the mission [km/s]

% Initialize times:
t_d = x(1); t_f = x(2); t_a = x(3); 

% Planets id numbers:
id_d = 1;                              % Mercury
id_f = 4;                              % Mars
id_a = 40;                             % Asteroid

[kep_d,ksun] = uplanet(t_d, id_d);
[rr_d,vv_d] = kep2car(kep_d(1),kep_d(2),kep_d(3),kep_d(4),kep_d(5),kep_d(6),ksun);
[kep_f,ksun] = uplanet(t_f, id_f);
[rr_f,vv_f] = kep2car(kep_f(1),kep_f(2),kep_f(3),kep_f(4),kep_f(5),kep_f(6),ksun);
kep_a = ephAsteroids(t_a, id_a);
[rr_a,vv_a] = kep2car(kep_a(1),kep_a(2),kep_a(3),kep_a(4),kep_a(5),kep_a(6),ksun);

[~,~,~,err1,vt1_i,vt1_f,~,~] = lambertMR(rr_d,rr_f,(t_f-t_d)*24*3600,ksun,0,0,0,0);

[~,~,~,err2,vt2_i,vt2_f,~,~] = lambertMR(rr_f,rr_a,(t_a-t_f)*24*3600,ksun,0,0,0,0);

dv = NaN; % Initialize output


% Check: there must be no errors in lambertMR
if (err1==0)&&(err2==0)
    
    vinf_m = vt1_f' - vv_f;
    vinf_p = vt2_i' - vv_f;
    
    if ~flag_atm 
    dvp  = PowerGravityAssist(vinf_m, vinf_p...
                    ,data.Mars.Radius, 0, data.Mars.mu, 1);
    elseif flag_atm
    dvp  = PowerGravityAssist(vinf_m, vinf_p...
        ,data.Mars.Radius, data.Mars.h_atm, data.Mars.mu, 1);
        if isnan(dvp)
            disp('!!! Warning: hyperbloic transfer not feasible')
        end
    end
    
    dv = norm(vv_d - vt1_i')+norm(vv_a - vt2_f') + dvp;
    
end



end
