function [logdet]=F_calc_det_QR(p, H, U)

[~,r]=size(U);
C = H*U;
if p <= r
    epsilon=1*10^(-8);
    logdet =det(C'*C+epsilon*eye(r,r));
else
    logdet = det(C'*C);
end
end
