function [Det_calculation]=calculation_det(p,r,H,Psi)
C=H*Psi;
if p <= r
    Det_calculation=log(det(C*C'));
else
    Det_calculation=log(det(C'*C));
end
end