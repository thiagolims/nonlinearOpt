% Make sure the workspace is clean before we start
clc
clear
clear global

addpath(genpath('library'));
addpath(genpath('algorithms'));
addpath(genpath('testFunctions'));

%% 
% intial point
x0=[-1; 2];

% Ghraphic: plot the function to be minimized
xmin = -2;
xmax = 2;
dx = .05;
ymin = -1;  
ymax = 3;
dy = .05;

[z,y]= meshgrid(xmin:dx:xmax, ymin:dy:ymax);
banana=10*z.^4 - 20*(z.^2).*y + 10*y.^2 + z.^2 -2*z + 5;
contour(z,y,banana,30);
axis('square');
hold;

%% Creating an automatic diff object
[x,y] = initVariablesADI(x0(1),x0(2));
xk = [x;y]; %% initial point with ADI

% f = fban(xk);
% g = gban(xk);


%% Anonymous functions 
fbanana = @(xk) fban(xk);
gbanana = @(xk) gban(xk);

fsquare = @(xk) squareX(xk);

%% Line Search Type
ls = 1; % Armijo = 1, Polynomial = 2

%% Invoking optimization algorithm
[xs, fs, gs] = gradientDescent(fbanana, gbanana, xk, ls);

