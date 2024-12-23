function [c, ceq] = nonlcon(x, data)
%
% nonlcon
% 
%
% PROTOTYPE:
%  dv = dvFun(x)
% 
% INPUT:
%  x [3,1]        Times array:         x(1) = departure time [MJD2000]
%                                      x(2) = first transfer tof [days]
%                                      x(3) = second transfer tof [days]
% 
% OUTPUT:
%  c              Array of nonlinear inequality ( c < 0 )
%  ceq            Array of nonlinear equality ( ceq = 0 )



% Initialize times:
x1 = cumsum(x);
t_d = x1(1); t_f = x1(2); t_a = x1(3); 

% Planets id numbers:
id_d = 1;                              % Mercury
id_f = 4;                              % Mars
id_a = 40;                             % Asteroid

[kep_d,ksun] = uplanet(t_d, id_d);
[rr_d,~] = kep2car(kep_d(1),kep_d(2),kep_d(3),kep_d(4),kep_d(5),kep_d(6),ksun);
[kep_f,ksun] = uplanet(t_f, id_f);
[rr_f,vv_f] = kep2car(kep_f(1),kep_f(2),kep_f(3),kep_f(4),kep_f(5),kep_f(6),ksun);
kep_a = ephAsteroids(t_a, id_a);
[rr_a,~] = kep2car(kep_a(1),kep_a(2),kep_a(3),kep_a(4),kep_a(5),kep_a(6),ksun);

[~,~,~,err1,~,vt1_f,~,~] = lambertMR(rr_d,rr_f,(t_f-t_d)*24*3600,ksun,0,0,0,0);

[~,~,~,err2,vt2_i,~,~,~] = lambertMR(rr_f,rr_a,(t_a-t_f)*24*3600,ksun,0,0,0,0);

rp = NaN; % Initialize output

% Check: there must be no errors in lambertMR
if (err1==0)&&(err2==0)
    
    vinf_m = vt1_f' - vv_f;
    vinf_p = vt2_i' - vv_f;
    
    [~,~,rp]  = PowerGravityAssist(vinf_m, vinf_p...
                    ,data.R_mars, 0, data.mu_mars);
    
    
end

ceq = 0;
tolerance = 1e-3;
c = data.R_mars + data.h - rp + tolerance;