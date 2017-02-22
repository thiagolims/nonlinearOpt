function [ gf ] = gban( x )
%GBAN Gradient of banana function
gf = x; %% init jacobian
gf(1)= 40.*x(1).^3-40.*x(1).*x(2)+2.*x(1)-2;
gf(2)= -20.*x(1).^2 + 20.*x(2);
gf = gf';
end

