function video_maker(num_video, Y, Itest, meansst, Psi, Xestimate, mask, time, r, p, filename,sensors)
v = VideoWriter(filename);
open(v);

     for i=1:num_video%length(time) 
             x = Y(:,Itest(i))-meansst;
        %     bounds = [min(x+meansst) max(x+meansst)];
              figure(1);           
            if p==1
               display_sensors_color_for_movie(x+meansst,mask,[], filename,time(i),r,p);
            else
                xls=Psi(:,1:r)*Xestimate(:,i);
                display_sensors_color_for_movie(xls+meansst,mask,sensors, filename,time(i),r,p);
            end
            frame = getframe(figure(1));
            writeVideo(v,frame);
             close(1)
     end
      close(v);
end 