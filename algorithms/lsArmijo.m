% -------------------------------------------------------------------
% METODO: LSArmijo
% Armijo Criterion for Line Search
% Author: Thiago Lima Silva (thiagolims@gmail.com)
% -------------------------------------------------------------------

function [ s, xs, fxs] = lsArmijo(f, x, d, g)
%LSARMIJO Armijo criteria for line search
c1 = 0.1;                   
alpha = 0.7;              % coefficient for backtracking

g0 = g;                   % f'(0)
f0 = f(x);                % f(0)
s  = 1;                   % initial step length

xs = x + s*d;             % x(s)
fxs = f(xs);              % f(xs)

gfd = g0'*d;              
if gfd > 0    
    xs = x;
    fxs = f0;
    s = 0;
    warning('No descent direction.')
    return;
end
    
while (fxs >  (f0 + s*c1*g0'*d))       %% Armijo condition
    s = alpha*s;        % new alpha
    xs = x + s*d;       % x(alpha)
    fxs = f(xs);    
end

