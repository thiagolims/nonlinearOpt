% -------------------------------------------------------------------
% METODO: gradient-based nonlinear optimization
% Newton, Quasi-Newton, Steepest Descent
% Author: Thiago Lima Silva (thiagolims@gmail.com)
% -------------------------------------------------------------------
function [ xs, fs, gs] = gradientDescent(f, g, h, x0, optType, lsType, varargin)
%gradientDescent steepest gradient descent algorithm
% input:
%   - f: pointer to the function to be minimized
%   - g: pointer to the function gradient
%   - h: point to the function hessian
%   - x0: initial point
%   - optType: direction search method (1 - Steepest Descent, 2 - Newton, 3 - DFP, 4 - BFGS, 5 - SR1/BFGS )
%   - lsType: line search type (1 - Armijo, 2 - ArmijoGoldstein, 3 - Polynomial)
% output:
%   - xs: local optimal solution
%   - fs: optimal function value f(xs)
%   - gs: gradient of fs at the optimal value df(xs)/dx

%% Initialization
xs = x0;
dim = size(xs,1);
f0 = f(x0);
g0 = g(x0);
h0 = h(x0);

fs = f(xs);
gs = g(xs)';
hs = h(xs);

d = -g0;
s = 1;

% BFGS and DFP parameters
H = eye(dim,dim);

% Parameter for Marquardt's Modification to the Newton Method
lambda = 1000;

%% convergence parameters
eps = 1.e-4;
tol = 1.e-6;
maxIter = 500;

if (norm(g0) < eps) return; end

for i=1:maxIter
    switch(optType)
        case 1, % Steepest Descent
            d = -gs;
            
            switch(lsType)
                case 1,
                    [s, x1, f1] = lsArmijo(f, double(xs), double(d), double(gs));
                    if s==0
                        error('zero step size determined.')
                    end
                    g1 = g(x1)';
                    %%TODO: compute gradient with ADI
                    
                case 2,
                    [s, x1, f1] = lsArmijoGoldstein(f, double(xs), double(d), double(gs));
                    g1 = g(x1)';
                    
                    
                case 3,
                    [s,x1, f1] = lsPolynomial(f, xs, d);
                    g1 = g(x1)';
            end
            
        case 2, % Newton
            % Parameters for Marquardt's Modification to the Newton Method
            beta = 5.0e-1*norm(hs)^-1;
            teta = 1e-6;
            
            if rcond(H)
                d = -(hs\gs);  % Newton's method direction
            else
                d = -(hs + lambda.*eye(dim,dim))\gs;   % Marquardt Modification to Newton's Method
            end
            
            while (gs'*d >= 0) ||  (gs'*d > -teta*norm(gs)*norm(d)) || (norm(d) < beta*norm(gs))
                lambda = 2*lambda;
                if rcond(hs + lambda.*eye(dim,dim)) > 1.0e-6
                    d = -(hs + lambda.*eye(dim,dim))\gs;
                else
                    d = -gs;
                end
            end
            
            switch(lsType)
                case 1,
                    [s, x1, f1] = lsArmijo(f, double(xs), double(d), double(gs));
                    if s==0
                        error('zero step size determined.')
                    end
                    g1 = g(x1)';
                    h1 = h(x1);
                    %%TODO: compute gradient with ADI
                    
                case 2,
                    [s, x1, f1] = lsArmijoGoldstein(f, double(xs), double(d), double(gs));
                    g1 = g(x1)';
                    h1 = h(x1);
                    
                    
                case 3,
                    [s,x1, f1] = lsPolynomial(f, xs, d);
                    g1 = g(x1)';
                    h1 = h(x1);
            end
            
            %% parameters update
            lambda = 0.5*lambda;
            
        case 3,  % DFP (Quasi-Newton)
            d = -H*gs; %direction
            
            switch(lsType)
                case 1,
                    [s, x1, f1] = lsArmijo(f, double(xs), double(d), double(gs));
                    if s==0
                        error('zero step size determined.')
                    end
                    g1 = g(x1)';
                    %%TODO: compute gradient with ADI
                    
                case 2,
                    [s, x1, f1] = lsArmijoGoldstein(f, double(xs), double(d), double(gs));
                    g1 = g(x1)';
                    
                    
                case 3,
                    [s,x1, f1] = lsPolynomial(f, xs, d);
                    g1 = g(x1)';
            end
            %% DFP Updating Formula
            y = g1-gs;
            sk = s*d; % x_{k+1} - x_{k}
            
            H = H - ((H*y)*(y'*H))/(y'*H*y) + (sk*sk')/(y'*sk);  %%  Nocedal book (page 139, eq 6.15)
            
        case 4,  % BFGS (Quasi-Newton)
            d = -H*gs; % direction
            
            switch(lsType)
                case 1,
                    [s, x1, f1] = lsArmijo(f, double(xs), double(d), double(gs));
                    if s==0
                        error('zero step size determined.')
                    end
                    g1 = g(x1)';
                    
                case 2,
                    [s, x1, f1] = lsArmijoGoldstein(f, double(xs), double(d), double(gs));
                    g1 = g(x1)';
                    
                    
                case 3,
                    [s,x1, f1] = lsPolynomial(f, xs, d);
                    g1 = g(x1)';
            end
            %% BFGS Updating Formula
            y = g1-gs;
            sk = s*d; %x_{k+1} - x_k
            
            H = bfgsUpdate(H, y, sk);            
        case 5,   %% SR1 (Quasi-Newton)
            dp = d;     % store previous direction
            sp = s;     % store previous step length
            
            d = -H*gs;  % next direction
            switch(lsType)
                %%Try SR1 update and switch to BFGS update in case a descent direction is not found.
                case 1,
                    [s, x1, f1] = lsArmijo(f, double(xs), double(d), double(gs));
                    if s==0
                        %% Step back to previous iterate before SR1 update
                        xp = xs - sp*dp;
                        gp = g(xp)';
                        H = Hp;
                        
                        y = gs - gp;
                        sk = sp*dp;
                                                
                        
                        %% BFGS Update
                        H = bfgsUpdate(H, y, sk);
                        
                        d = -H*gs;  % new direction base on a rank-2 Hessian approximation
                        
                        [s, x1, f1] = lsArmijo(f, double(xs), double(d), double(gs));
                    end
                    g1 = g(x1)';
                    
                case 2,
                    [s, x1, f1] = lsArmijoGoldstein(f, double(xs), double(d), double(gs));
                    if s==0
                        %% Step back to previous iterate before SR1 update
                        xp = xs - sp*dp;
                        gp = g(xp)';
                        H = Hp;
                        
                        y = gs - gp;
                        sk = sp*dp;
                                                
                        
                        %% BFGS Update
                        H = bfgsUpdate(H, y, sk);
                        
                        d = -H*gs;  % new direction base on a rank-2 Hessian approximation
                        [s, x1, f1] = lsArmijoGoldstein(f, double(xs), double(d), double(gs));
                    end
                    g1 = g(x1)';
                    
                    
                case 3,
                    [s,x1, f1] = lsPolynomial(f, xs, d);
                    g1 = g(x1)';
            end
            
            y = g1-gs; % df/dx_{k+1} - df/dx_k
            sk = s*d;  % x_{k+1} - x_k
            
            
            % if denR1 is small, keep the same hessian approximation H
            denR1 = ((sk - H*y)'*y);
            r = 1e-08;
            if (sign(denR1)*denR1) >= (r*norm(y)*norm(sk - H*y))
                Hp = H;  %% store previous Hessian approximation
                H = H + (sk - H*y)*(sk - H*y)'/denR1; % unique rank-one updating formula satisfying the secant equation
            end
            
    end
    
    
    % Convergence test
    gradientNorm = (norm(double(g1))/norm(double(g0)));
    
    
    if (gradientNorm < eps)
        disp(sprintf('\n Iteration number %d', i));
        disp(sprintf('  f = %g  Step = %d', fs, s));
        disp(sprintf('  X = ')); disp(sprintf(' %d ',xs));
        disp(sprintf('  g/g0 = ')); disp(sprintf(' %d ',gradientNorm));
        
        disp('Convergence OK!');
        return
    end
    
    % Plotting and printing
    xsD = double(xs);
    x1D = double(x1);
    
    dx = [xsD(1) x1D(1)];
    dy = [xsD(2) x1D(2)];
    plot(dx,dy,'-ro');
    
    disp(sprintf('\n Iteration number %d', i));
    disp(sprintf('  f = %g  Step = %d', f1, s));
    disp(sprintf('  X = ')); disp(sprintf(' %d ',x1));
    disp(sprintf('  g/g0 = ')); disp(sprintf(' %d ',gradientNorm));
    
    xs = x1;
    gs = g1;
    fs = f1;
    if optType == 2 % Newton
        hs = h1;
    end
    
end

end

