function [Perror_ave,Perror_std]=make_ensopod(r,Xorg,Psi,S,V,x,num_video, Y, Itest, meansst, mask, time,TT)
p=1;
filename=('enso_POD');
Xestimate = (Psi(:,1:r)'*x);
filename=['enso_proj_mode',num2str(r)];
filename = ['movie', '\',filename, '.avi'];
video_maker(num_video, Y, Itest, meansst, Psi, Xestimate, mask, time, r, p, filename);

for ii=1:TT
    xes=Psi(:,1:r)*S(1:r,1:r)*V(ii,1:r)';
    Perror(ii) = norm(xes-Xorg(:,ii))/norm(Xorg(:,ii));
end
Perror_ave = mean(Perror);
Perror_std = std(Perror);
end