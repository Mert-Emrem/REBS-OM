function E = KeplerSolver(t, e, a, mu, t0, E0, tol)

n = sqrt(mu/a^3);
T = 2*pi/n;
k = (t-t0)/T;

f = @(E) E - e.*sin(E) - n.*t;

E_guess = @(t, e) n.*t + e.*sin(n.*t)./(1-sin(n.*t+e)+sin(n.*t));

if nargin < 7
    tol = 1e-6;
end 

options = optimoptions('fsolve', 'TolFun', tol, 'TolX', tol);
E = fsolve(f, E_guess(t,e)-E0, options);

E = k*2*pi + E;


end