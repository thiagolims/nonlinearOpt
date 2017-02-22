function [ xs, fs, gs] = gradientDescent(f, g, x0, ls)
%gradientDescent steepest gradient descent algorithm
% input:
%   - f: pointer to the function to be minimized
%   - g: pointer to the gradient of func
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
g0 = [f0.jac{1}; f0.jac{2}];

fs = f(xs);
gs = g(xs)';
% jacFs = [fs.jac{1}; fs.jac{2}];
% gs = fs;
% gs.val = jacFs;

%% convergence parameters
eps = 1.e-4; 
tol = 1.e-6;
maxIter = 1500;

%%TODO: implement the steepest descent algorithm
if (norm(g0) < eps) return; end

for i=1:maxIter
    d = -gs;
    switch(ls)
        case 1,
            [s, x1, f1] = lsArmijo(f, double(xs), double(d), double(gs));            
            g1 = g(x1)';
            %%TODO: compute gradient with ADI
            
%             x1 = xs + s*d;
%             f1 = f(x1);
%             g1 = g(x1);
%             jacF1 = [f1.jac{1}; f1.jac{2}];
%             g1 = f1;
%             g1.val = jacF1;           
         
        case 2,
            s = 1; %% full step  %%TODO: implement polynomial
            x1 = xs + s*d;
            f1 = f(x1);
            g1 = g(x1);
%             jacF1 = [f1.jac{1}; f1.jac{2}];
%             g1 = f1;
%             g1.val = jacF1;    
           
    end
    
     % Convergence test
     gradientNorm = (norm(double(g1))/norm(double(gs)));
     
     if (gradientNorm < eps)
         disp(sprintf('\n Iteration number %d', i));
         disp(sprintf('  f = %g  Step = %d', f, s));
         disp(sprintf('  X = ')); disp(sprintf(' %d ',xs));
         disp(sprintf('  g/g0 = ')); disp(sprintf(' %d ',gradientNorm));   

     	disp('Convergence OK!');
        return 
     end
     
    % Plotting and printing
    xsD = double(xs);
    x1D = double(xs);
    
    dx = [xsD(1) x1D(1)];
    dy = [xsD(2) x1D(2)];
    plot(dx,dy,'-ro');
    
%     disp(sprintf('\n Iteration number %d', i));
% 	disp(sprintf('  f = %g  Step = %d', f1, s));
% 	disp(sprintf('  X = ')); disp(sprintf(' %d ',X1));
% 	disp(sprintf('  g/g0 = ')); disp(sprintf(' %d ',norma));   
    xs = x1;
    gs = g1;
    fs = f1;
end

end

