function R = iner2perifocal(om, i, OM)

R_om = [cos(om), sin(om), 0; -sin(om), cos(om), 0; 0, 0, 1];
R_i = [1, 0, 0; 0, cos(i), sin(i); 0, -sin(i), cos(i)];
R_OM = [cos(OM), sin(OM), 0; -sin(OM), cos(OM), 0; 0, 0, 1];

R = R_om * R_i * R_OM;

end