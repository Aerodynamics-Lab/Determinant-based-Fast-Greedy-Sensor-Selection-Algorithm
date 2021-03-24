function makevideo_NOAASST(r, num_video, x, Y, Itest, meansst, Psi, mask, time, p, ...
    sensors_rand, sensors_convex, sensors_QR, sensors_DG, sensors_QD, ...
    Xestimate_rand, Xestimate_convex, Xestimate_QR,Xestimate_DG,Xestimate_QD)

%% Rand
name='rand';
filename=['enso_', name ,'_mode' num2str(r) '_sensor',num2str(p)];
display_sensors_color(x(:,1)+meansst,mask,sensors_rand, filename);

filename = ['movie', '\',filename, '.avi'];
video_maker(num_video, Y, Itest, meansst, Psi, Xestimate_rand, mask, time, r, p, filename, sensors_rand);  
%% convex
 name='convex';
filename=['enso_', name ,'_mode' num2str(r) '_sensor',num2str(p)];
display_sensors_color(x(:,1)+meansst,mask,sensors_convex, filename);  
filename = ['movie', '\',filename, '.avi'];
video_maker(num_video, Y, Itest, meansst, Psi, Xestimate_convex, mask, time, r, p, filename, sensors_convex);
%% QR
 name='QR';
filename=['enso_', name ,'_mode' num2str(r) '_sensor',num2str(p)];
display_sensors_color(x(:,1)+meansst,mask,sensors_QR, filename);  
filename = ['movie', '\',filename, '.avi'];
video_maker(num_video, Y, Itest, meansst, Psi, Xestimate_QR, mask, time, r, p, filename, sensors_QR);
%% DG
 name='DG';
filename=['enso_', name ,'_mode' num2str(r) '_sensor',num2str(p)];
display_sensors_color(x(:,1)+meansst,mask,sensors_DG, filename);
filename = ['movie', '\',filename, '.avi'];
video_maker(num_video, Y, Itest, meansst, Psi, Xestimate_DG, mask, time, r, p, filename, sensors_DG);
%% QD
 name='QD';
filename=['enso_', name ,'_mode' num2str(r) '_sensor',num2str(p)];
display_sensors_color(x(:,1)+meansst,mask,sensors_QD, filename);
filename = ['movie', '\',filename, '.avi'];
video_maker(num_video, Y, Itest, meansst, Psi, Xestimate_QD, mask, time, r, p, filename, sensors_QD);
end