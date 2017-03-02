function [H1] = bfgsUpdate(H, y, sk)
%bfgsUpdate BFGS updating formula
pk = 1/(y'*sk);  %% Nocedal Book (Page 139, eq. 6.14)   
H1 = (eye(length(sk)) - pk*sk*y')*H*(eye(length(sk)) - pk*y*sk') + pk*(sk*sk'); %% Nocedal Book (Page 140, eq. 6.17)
end

