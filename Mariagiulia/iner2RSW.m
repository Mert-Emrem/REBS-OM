function R = iner2RSW(om, th, i, OM)

R_omth = [cos(om+th), sin(om+th), 0; -sin(om+th), cos(om+th), 0; 0, 0, 1];
R_i = [1, 0, 0; 0, cos(i), sin(i); 0, -sin(i), cos(i)];
R_OM = [cos(OM), sin(OM), 0; -sin(OM), cos(OM), 0; 0, 0, 1];

R = R_omth * R_i * R_OM;

end