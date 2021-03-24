function [H]=makesensor_matrix(p, n, isensors)
H=zeros(p,n);    %Create array of all zeros
    for t=1:p
       H(t, isensors(t))=1;
    end