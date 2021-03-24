function [time_QR, H, sensors]=sensor_QR(Psi, r, p, n)
tic;
[sensors]=subsensor_QR(r, p, Psi(:,1:r));
time_QR=toc;
[H]=makesensor_matrix(p, n, sensors);
end