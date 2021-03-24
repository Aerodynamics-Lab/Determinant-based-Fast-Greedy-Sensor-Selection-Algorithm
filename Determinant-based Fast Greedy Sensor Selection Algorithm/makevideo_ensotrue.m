function makevideo_ensotrue(num_video, Y, Itest, meansst, Psi,  mask, time, r)
%% select validation snapshot
 p=1;
 filename=('enso_true');
 filename = ['movie', '\',filename, '.avi'];
video_maker(num_video, Y, Itest, meansst, Psi, [], mask, time, r, p, filename);
end