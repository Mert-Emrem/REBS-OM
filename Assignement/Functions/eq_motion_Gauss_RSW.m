function ds = eq_motion_Gauss_RSW(t, kep_el, acc_pert_function,  mu) %, J2, Re, day)

a = kep_el(1);
e = kep_el(2); 
i = kep_el(3); 
OM = kep_el(4); 
om = kep_el(5);
th = kep_el(6);


%a_p_RSW = acc_pert_function_RSW(t, kep_el); % ) 
a_p = acc_pert_function(t, kep_el);  

R = iner2RSW(om,th, i, OM);
a_p_RSW = R*a_p;
a_r = a_p_RSW(1);
a_s = a_p_RSW(2);
a_w = a_p_RSW(3);

p = a*(1-e^2);
h = sqrt(p*mu);
r = p/(1+e*cos(th));

da = 2*a^2/h*(e*sin(th)*a_r + p/r * a_s);
de = 1/h * (p*sin(th)*a_r + ((p+r)*cos(th) + r*e)*a_s);
di = r*cos(th+om)/h * a_w;
dOM = r*sin(th+om)/(h*sin(i)) * a_w;
dom = 1/(h*e) * (-p*cos(th)*a_r + (p+r)*sin(th)*a_s) - r*sin(th+om)*cos(i)/(h*sin(i))*a_w;
dth = h/r^2 + 1/(e*h)*(p*cos(th)*a_r - (p+r)*sin(th)*a_s);

ds = [da; de; di; dOM; dom; dth];

end


