function [ f ] = fban(x)
%FBAN Banana function
f=10*x(1).^4 - 20*x(1).^2*x(2) + 10*x(2).^2 + x(1).^2 -2*x(1) + 5;
end

