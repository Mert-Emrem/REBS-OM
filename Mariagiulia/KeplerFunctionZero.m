function E = KeplerFunctionZero(t, e, a, mu, t0, E0, tol)

if nargin <7
    tol = 1e-6;
end

n = sqrt(mu/a^3);
T = 2*pi/n;
k = floor((t-t0)/T);
t_bar = t - k*T;


f = @(E) E - e.*sin(E) - E0 + e.*sin(E0) - n.*(t_bar - t0);
E_guess = @(t, e) n.*t + (e.*sin(n.*t))./(1-sin(n.*t + e) + sin(n.*t));

%options = optimoptions('fsolve', 'TolFun', tol, 'TolX', tol);

E = fzero(f, E_guess(t_bar, e));

E = E + k*2*pi;
