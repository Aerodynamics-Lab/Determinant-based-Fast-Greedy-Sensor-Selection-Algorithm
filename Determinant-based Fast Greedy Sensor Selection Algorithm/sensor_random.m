function [time_rand, H,sensors]=sensor_random(n,p)
tic;
sensors = randperm(n,p);
time_rand=toc;

[H]=makesensor_matrix( p, n, sensors);
end