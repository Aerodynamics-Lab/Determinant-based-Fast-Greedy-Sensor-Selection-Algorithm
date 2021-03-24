function [Normalized_det,det,time]=data_organization(q,v,time_rand,time_convex,time_QR,time_DG,time_QD,det_rand,det_convex,det_QR,det_DG,det_QD)

det(1:v,1)= 1:q;
det(1:v,2)=det_rand(1:v,1);
det(1:v,3)=det_convex(1:v,1);
det(1:v,4)=det_QR(1:v,1);
det(1:v,5)=det_DG(1:v,1);
det(1:v,6)=det_QD(1:v,1);

time(1:v,1)= 1:q;
time(1:v,2)=time_rand(1:v,1);
time(1:v,3)=time_convex(1:v,1);
time(1:v,4)=time_QR(1:v,1);
time(1:v,5)=time_DG(1:v,1);
time(1:v,6)=time_QD(1:v,1);

for z=1:v
    Normalized_det(z,1)=z;
    Normalized_det(z,2)=det_rand(z,1)/det_QR(z,1);
    Normalized_det(z,3)=det_convex(z,1)/det_QR(z,1);
    Normalized_det(z,4)=det_QR(z,1)/det_QR(z,1);
    Normalized_det(z,5)=det_DG(z,1)/det_QR(z,1);
    Normalized_det(z,6)=det_QD(z,1)/det_QR(z,1);
end

end