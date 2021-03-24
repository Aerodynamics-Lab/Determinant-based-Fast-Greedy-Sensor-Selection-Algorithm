clear; close all;
warning('off','all')

%% Selection of Problems !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%num_problem=1; %%Randomized sensor matrix
num_problem=2; %%NOAA-SST
% ^^^ It takes a long time to obtain the solution in the convex approximation method.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%% Initinal Values
r=10;
num=r+10;
maxiteration=200; % for convex
n=2000; %for Randomized sensor matri
num_ave=1;
TT=52*10;%10-years;
num_video=1;% maxmum:TT

if  num_problem == 2%%NOAA-SST
    [Psi,S,V,Xorg,Itest,mask, meansst, Y, x, time]=read_NOAA_SST(TT); %Reading a NOAA-SST
    [n,~]=size(Psi);
    %Need movie?
    mkdir('movie');
    makevideo_ensotrue(num_video, Y, Itest, meansst, Psi,  mask, time, r) % select validation snapshot
end

v=0;
for q=1:1%num %Parameter loop
    v=v+1; %Counter
    p=q; %Parameter selection
    
    if  num_problem == 2 %NOAA-SST
        if v == 1
            [Perror_ave_pod,Perror_std_pod]=make_ensopod(r,Xorg,Psi,S,V,x,num_video, Y, Itest, meansst, mask, time,TT);% POD approximation with r eigenssts
        end
    end
    
    for w=1:num_ave   % Average loop
        
        if num_problem == 1 %Randomized sensor matrix
            Psi=randn(n,r);
        end
        
        %% Sparse sensor matrix
        % Random selection
        [time_rand(v,w+1), H_rand,sensors_rand]=sensor_random(n,p);
        [det_rand(v,w+1)]=calculation_det(p,r,H_rand,Psi);
        % convex approximation
        %!! I recommend you use the following dummy values if you do not need the solution in the convex approximation in NOAA-SST.
        %  [time_convex(v,w+1), H_convex,sensors_convex]=sensor_convex(Psi, r, p, n, maxiteration);
        %   [det_convex(v,w+1)]=calculation_det(p,r,H_convex,Psi);
        %***************************************
        time_convex(v,w+1)=time_rand(v,w+1);
        det_convex(v,w+1)=det_rand(v,w+1);
        H_convex=H_rand;
        sensors_convex=sensors_rand;
        %***************************************
        % QR
        [time_QR(v,w+1), H_QR, sensors_QR]=sensor_QR(Psi, r, p, n);
        [det_QR(v,w+1)]=calculation_det(p,r,H_QR,Psi);
        % DG
        [time_DG(v,w+1), H_DG, sensors_DG]=sensor_DG(Psi, r, p, n);
        [det_DG(v,w+1)]=calculation_det(p,r,H_DG,Psi);
        % QD
        [time_QD(v,w+1), H_QD, sensors_QD]=sensor_QD(Psi, r, p, n);
        [det_QD(v,w+1)]=calculation_det(p,r,H_QD,Psi);
        
        %% Percent error calculation
        if  num_problem == 2 %NOAA-SST
            [Xestimate_rand, Xestimate_convex, Xestimate_QR,Xestimate_DG,Xestimate_QD,...
                Perror_rand(v,w+1),Perror_std_rand(v,w+1),Perror_convex(v,1),Perror_std_convex(v,1),...
                Perror_QR(v,1),Perror_std_QR(v,1),Perror_DG(v,1),Perror_std_DG(v,1),Perror_QD(v,1),Perror_std_QD(v,1)]=...
                calculation_error(v, q, w, r, TT, Xorg, Psi, H_rand, H_convex, H_QR, H_DG, H_QD);
            Perror_rand(v,1)=mean(Perror_rand(v,2:w+1));
            Perror_std_rand(v,1)=mean(Perror_std_rand(v,2:w+1));
        end
        
    end
    sensor_memo=zeros(p,5);
    sensor_memo(1:p,1)=sensors_rand(1:p)';
    sensor_memo(1:p,2)=sensors_convex(1:p);
    sensor_memo(1:p,3)=sensors_QR(1:p)';
    sensor_memo(1:p,4)=sensors_DG(1:p);
    sensor_memo(1:p,5)=sensors_QD(1:p)';
    filename=['sensors_p_',num2str(p),'.mat'];
    save(filename,'sensor_memo');
    if  num_problem == 2 %NOAA-SST
        % Arrange
        Perror(v,1)= q;
        Perror(v,2)= Perror_rand(v,1);
        Perror(v,3)= Perror_std_rand(v,1);
        Perror(v,4)= Perror_convex(v,1);
        Perror(v,5)= Perror_std_convex(v,1);
        Perror(v,6)= Perror_QR(v,1);
        Perror(v,7)= Perror_std_QR(v,1);
        Perror(v,8)= Perror_DG(v,1);
        Perror(v,9)= Perror_std_DG(v,1);
        Perror(v,10)= Perror_QD(v,1);
        Perror(v,11)= Perror_std_QD(v,1);
        Perror(v,12)= Perror_ave_pod;
        % Save
        save('Perror_all.mat','Perror');
        save('Perror_rand.mat','Perror_rand');
        save('Perror_convex.mat','Perror_convex');
        save('Perror_QR.mat','Perror_QR');
        save('Perror_DG.mat','Perror_DG');
        save('Perror_QD.mat','Perror_QD');
    end
    %% Average operation&Save
    [time_rand,time_convex,time_QR,time_DG,time_QD,det_rand,det_convex,det_QR,det_DG,det_QD]...
        =averaged_operation(v,w,time_rand,time_convex,time_QR,time_DG,time_QD,det_rand,det_convex,det_QR,det_DG,det_QD);
    %% Video
    if  num_problem == 2 %NOAA-SST
        makevideo_NOAASST(r, num_video, x, Y, Itest, meansst, Psi, mask, time, p, ...
            sensors_rand, sensors_convex, sensors_QR, sensors_DG, sensors_QD, ...
            Xestimate_rand, Xestimate_convex, Xestimate_QR,Xestimate_DG,Xestimate_QD)
    end
    %% Display
    %Display(v,  w, time_rand, time_convex, time_QR, time_DG, time_QD, det_rand, det_convex, det_QR, det_DG, det_QD);
    [Normalized_det,det,time_all]...
        =data_organization(q,v,time_rand,time_convex,time_QR,time_DG,time_QD,det_rand,det_convex,det_QR,det_DG,det_QD);
    text=['CalTime_p:',num2str(p),'_r:',num2str(r),'_n:',num2str(n),'_RA_',num2str(time_all(v,2)),'_CO_',num2str(time_all(v,3)),...
        '_QR_',num2str(time_all(v,4)),'_DG_',num2str(time_all(v,5)),'_QD_',num2str(time_all(v,6))];
    disp(text);
    text=['Determinat_p:',num2str(p),'_r:',num2str(r),'_n:',num2str(n),'_RA_',num2str(Normalized_det(v,2)),'_CO_',num2str(Normalized_det(v,3)),...
        '_QR_',num2str(Normalized_det(v,4)),'_DG_',num2str(Normalized_det(v,5)),'_QD_',num2str(Normalized_det(v,6))];
    disp(text);
end
save('time_all.mat','time_all');
save('det_random.mat','det');
save('Normalized_det_random.mat','Normalized_det');

warning('on','all')
%         time_convex(v,w+1)=time_rand(v,w+1);
%         det_convex(v,w+1)=det_rand(v,w+1);
%         H_convex=H_rand;
%         sensors_convex=sensors_rand;
