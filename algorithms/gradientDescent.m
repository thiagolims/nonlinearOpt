% -------------------------------------------------------------------
% METODO: gradient-based nonlinear optimization
% Newton, Quasi-Newton, Steepest Descent
% Author: Thiago Lima Silva (thiagolims@gmail.com)
% -------------------------------------------------------------------
function [ xs, fs, gs] = gradientDescent(f, g, h, x0, ls)
%gradientDescent steepest gradient descent algorithm
% input:
%   - f: pointer to the function to be minimized
%   - g: pointer to the function gradient
%   - h: point to the function hessian
%   - x0: initial point
%   - ls: line search type (1 - Armijo, 2 - Polynomial)
% output:
%   - xs: local optimal solution
%   - fs: optimal function value f(xs)
%   - gs: gradient of fs at the optimal value df(xs)/dx

%% Initialization
xs = x0;
dim = size(xs,1);
H = eye(dim,dim);
f0 = f(x0);
g0 = g(x0);
h0 = h(x0);

fs = f(xs);
gs = g(xs)';
hs = h(xs);

%% convergence parameters
eps = 1.e-4; 
tol = 1.e-6;
maxIter = 500;

%%TODO: implement the steepest descent algorithm
if (norm(g0) < eps) return; end

for i=1:maxIter
%     d = -gs;
    d = -(hs\gs); % Newton
    switch(ls)
        case 1,
            [s, x1, f1] = lsArmijo(f, double(xs), double(d), double(gs));            
            g1 = g(x1)';
            h1 = h(x1);
            %%TODO: compute gradient with ADI
         
        case 2,
            s = 1; %% full step  %%TODO: implement polynomial
            x1 = xs + s*d;
            f1 = f(x1);
            g1 = g(x1);
           
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
    hs = h1;
end

end

