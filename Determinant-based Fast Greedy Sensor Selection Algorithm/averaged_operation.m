function [time_rand,time_convex,time_QR,time_DG,time_QD,...
    det_rand,det_convex,det_QR,det_DG,det_QD]...
    =averaged_operation(v,w,time_rand,time_convex,time_QR,time_DG,time_QD,...
    det_rand,det_convex,det_QR,det_DG,det_QD)

   %Random
   time_rand(v,1)=mean(time_rand(v,2:w+1));
   det_rand(v,1)=mean(det_rand(v,2:w+1));
   %Convex
   time_convex(v,1)=mean(time_convex(v,2:w+1));
   det_convex(v,1)=mean(det_convex(v,2:w+1));
   %QR
   time_QR(v,1)=mean(time_QR(v,2:w+1));
   det_QR(v,1)=mean(det_QR(v,2:w+1));
   %DG
   time_DG(v,1)=mean(time_DG(v,2:w+1));
   det_DG(v,1)=mean(det_DG(v,2:w+1));
   %QD
   time_QD(v,1)=mean(time_QD(v,2:w+1));
   det_QD(v,1)=mean(det_QD(v,2:w+1));
   
   %% Save
    save_operation(time_rand,time_convex,time_QR,time_DG,time_QD,det_rand,det_convex,det_QR,det_DG,det_QD);
end