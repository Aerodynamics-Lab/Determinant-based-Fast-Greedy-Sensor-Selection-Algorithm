function [time_QD, H,sensors]=sensor_QD(Psi, r, p, n)
if p <= r
    tic;
    [sensors]=subsensor_QR(r, p, Psi(:,1:r));
    time_QD=toc;
    [H]=makesensor_matrix(p, n, sensors);
else
    tic;
    [isensors]=subsensor_QR(r, r, Psi(:,1:r));
    [H]=makesensor_matrix(r, n, isensors);
    [sensors]=subsensor_DG_p(Psi(:,1:r), r, p, H, isensors');
    [H]=makesensor_matrix(p, n, sensors);
    time_QD=toc;
end
end