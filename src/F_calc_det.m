function [logdet]=F_calc_det(p, H, U)

    [~,r]=size(U);
    C = H*U;
    logdet = det(C'*C);
%     if p <= r
%         logdet = det(C*C');
%     else
%         logdet = det(C'*C);
%     end
end
