function [time_DG, H,sensors]=sensor_DG(Psi, r, p, n)

if p <= r
    tic;
    [sensors]= subsensor_DG_r(Psi(:,1:r), r, p);
    time_DG=toc;
    [H]=makesensor_matrix(p, n, sensors);
else
    tic;
    [isensors]= subsensor_DG_r(Psi(:,1:r), r, r);
    [H]=makesensor_matrix(r, n, isensors);
    [sensors]=subsensor_DG_p(Psi(:,1:r), r, p, H, isensors);
    time_DG=toc;
    [H]=makesensor_matrix(p, n, sensors);
end

end