function [rr, vv] = kep2car_vec(a,e,i,OM,om,th, mu)

% VUOLE VETTORI RIGA IN INPUT

% i = rad2deg(i);
% OM = rad2deg(OM);
% om = rad2deg(om);
% th = rad2deg(th);

p = a.*(1-e.^2);
h = sqrt(p.*mu);
r = p ./ (1+e.*cos(th));

r_PF = r .* [cos(th); sin(th); 0*sin(th)];
v_PF = (mu./h) .* [-sin(th); (e+cos(th)); 0*sin(th)];


% Rotation matrices: Earth-Centered Inertial --> Perifocal   (ECI->PF)

R_om = zeros(3, 3, length(om));

c_om = cos(om);
s_om = sin(om);

R_om(1, 1, :) = c_om;
R_om(1, 2, :) = s_om;
R_om(2, 1, :) = -s_om;
R_om(2, 2, :) = c_om;
R_om(3, 3, :) = 1;

R_i = zeros(3, 3, length(i));

c_i = cos(i);
s_i = sin(i);

R_i(1, 1, :) = 1;
R_i(2, 2, :) = c_i;
R_i(2, 3, :) = s_i;
R_i(3, 2, :) = -s_i;
R_i(3, 3, :) = c_i;

R_OM = zeros(3, 3, length(OM));

c_OM = cos(OM);
s_OM = sin(OM);
       
R_OM(1, 1, :) = c_OM;
R_OM(1, 2, :) = s_OM;
R_OM(2, 1, :) = -s_OM;
R_OM(2, 2, :) = c_OM;
R_OM(3, 3, :) = 1;
    
R313 = pagemtimes(R_om, R_i); % ECI --> PF
R313 = pagemtimes(R313, R_OM);

R313_t = pagetranspose(R313);

% l = R313(1,2, :);
% R313(1,2, :) = R313(2,1,:);
% R313(2,1,:) = l;
% 
% l = R313(1,3, :);
% R313(1,3, :) = R313(3,1,:);
% R313(3,1,:) = l;
% 
% l = R313(2,3, :);
% R313(2,3, :) = R313(3,2,:);
% R313(3,2,:) = l;
% 
% R313_t = R313;


% PF --> ECI
r_ECI = pagemtimes(R313_t, r_PF);
rr = r_ECI(:,:,1);
rr = squeeze(rr);
v_ECI = pagemtimes(R313_t, v_PF);
vv = v_ECI(:,:,1);
vv = squeeze(vv);


end
