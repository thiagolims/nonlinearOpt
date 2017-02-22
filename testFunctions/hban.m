function [ hf ] = hban( x )
%HBAN hessian of banana function
hf(1,1)= 120.*x(1).^2 - 40.*x(2) + 2; 
hf(2,2)= 20;

hf(1,2)= -40.*x(1); 
hf(2,1)= -40.*x(1);
end

