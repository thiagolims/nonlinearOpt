function [ s, xs, fxs] = lsArmijo(f, x, d, g)
%LSARMIJO Armijo criteria for line search

rho1 = 0.8;
rho2 = 0.1;
nUP = 2;                  % coefficient for interval increment
nDOWN = 0.7;              % coefficient for backtracking


g0 = g;                   % f'(0)
f0 = f(x);                % f(0)
s  = 1;                   % initial step length

xs = x + s*d;             % x(s)
fxs = f(xs);              % f(xs)

gfd = g0'*d;              
if gfd > 0
    error('No descent direction')
end
    
while (fxs <= (f0 + s*rho1*g0'*d)  ) %% Armijo condition
    s = nUP*s;
    xs = x + s*d;
    fxs = f(xs);
end

while (fxs >  (f0 + s*rho2*g0'*d))
    s = nDOWN*s;        % new alpha
    xs = x + s*d;       % x(alpha)
    fxs = f(xs);    
end

