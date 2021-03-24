function [isensors,zhat]=subsensor_convex(U, r, p, maxiteration)
[n,~]=size(U);
UU=[];

UU(:,:,1)=U(1:n,1:r); 

[zhat, ~, ~] = sens_sel_approxnt_vec(UU, p, maxiteration);
isensors = find(zhat);