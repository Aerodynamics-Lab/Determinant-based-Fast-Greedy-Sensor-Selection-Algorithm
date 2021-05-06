%% Main program
%% ///////////////////////////////////////////////////////////////////
% Comments:
% 	Collaborator: Yuji Saito, Keigo Yamada, Taku Nonomura
%                 Kumi Nakai, Takayuki Nagata
% 	Last modified: 2020/7/17
% Nomenclature:
% - Scalars
%   n : Number of degrees of freedom of spatial POD modes (state dimension)
%   p : Number of sensors
%   r : Number of rank for truncated POD
%   m : Number of snaphot (temporal dimension)
% - Matrices
% 	X : Supervising data matrix
% 	Y : Observation matrix
% 	H : Sparse sensor location matrix
% 	U : Spatial POD modes matrix
% 	C : Measurement matrix
% 	Z : POD mode amplitude matrix
%% ===================================================================

clear; close all;
warning('off','all')

%% Selection of Problems ============================================
% num_problem=1; % //Randomized sensor problem//
num_problem=2; % //NOAA-SST//
% !<NOAA-SST> It takes a long time to obtain the solution in the convex
% !<NOAA-SST> approximation method and the convex method is commented out
% !<NOAA_SST> as default setting for reduction of demo time.
%
%% Parameters =======================================================
r = 10;
pmin = 1;
pinc = 1;
pmax = 10;
ps   = pmin:pinc:pmax;
%ps = [10 15 20];
%ps = [1 2 3 4 5 6 7 8 9 10 12 14 16 18 20];
% ps = [1 2 3 4 5 6 7 8 9 10];
num_ave = 1; % Number of iteration for averaging operation
CNT = 0; % Counter
CV=5;
maxiteration = 200; % Max iteration for convex approximation
% //Randomized sensor problem//
n = 2000;
% //NOAA-SST//
m = 52*10; % 10-years (52weeks/year)
num_video = 1; % maxmum: m

%% Preparation of output directories ================================
workdir   = ('../work');
videodir  = [workdir,'/video'];
sensordir = [workdir,'/sensor_location'];
mkdir(workdir);
mkdir(videodir);
mkdir(sensordir);

%% Randomized sensor problem ========================================
if num_problem == 1
    
    %% Sensor selection =============================================
    for p = ps
        CNT = CNT+1;
        text = [ num2str(p),' sensor selection started --->' ];
        disp(text);
        
        %% Average loop =============================================
        for w=1:1:num_ave
            
            %% Preprocess for Randomized problem ====================
            U = randn(n,r);
            
            %% Random selection -------------------------------------
            [time_rand(CNT,w+1), H_rand, sensors_rand] = F_sensor_random(n,p);
            det_rand (CNT,w+1) = F_calc_det  (p,H_rand,U);
            tr_rand  (CNT,w+1) = F_calc_trace(p,H_rand,U);
            eig_rand (CNT,w+1) = F_calc_eigen(p,H_rand,U);
            
            %% D-optimality - Convex---------------------------------
            %!! This is very time consuming proceduce, We do not recommend to try this
            % [time_DC(CNT,w+1), H_DC, sensors_DC, ...
            %  NT_TOL_cal_DC(CNT,w+1), iter_DC(CNT,w+1)] ...
            %  = F_sensor_DC(U,p,maxiteration);
            % det_DC (CNT,w+1) = F_calc_det  (p,H_DC,U);
            % tr_DC  (CNT,w+1) = F_calc_trace(p,H_DC,U);
            % eig_DC (CNT,w+1) = F_calc_eigen(p,H_DC,U);
            %!! I recommend you use the following dummy values
            %   if you do not need the solution in the convex approximation in NOAA-SST.
            time_DC(CNT,w+1) = time_rand(CNT,w+1);
            det_DC (CNT,w+1) = det_rand (CNT,w+1);
            tr_DC  (CNT,w+1) = tr_rand  (CNT,w+1);
            eig_DC (CNT,w+1) = eig_rand (CNT,w+1);
            H_DC=H_rand;
            sensors_DC=sensors_rand;
            NT_TOL_cal_DC(CNT,w+1)=0;
            iter_DC(CNT,w+1)=0;
            
            %% Maximization of row norm - Greedy based on QR --------
            [time_QR(CNT,w+1), H_QR, sensors_QR] = F_sensor_QR(U,p);
            det_QR (CNT,w+1) = F_calc_det  (p,H_QR,U);
            tr_QR  (CNT,w+1) = F_calc_trace(p,H_QR,U);
            eig_QR (CNT,w+1) = F_calc_eigen(p,H_QR,U);
            
            %% D-optimality - Greedy --------------------------------
            [time_DG(CNT,w+1), H_DG, sensors_DG] = F_sensor_DG(U,p);
            det_DG (CNT,w+1) = F_calc_det  (p,H_DG,U);
            tr_DG  (CNT,w+1) = F_calc_trace(p,H_DG,U);
            eig_DG (CNT,w+1) = F_calc_eigen(p,H_DG,U);
            
            %% D-optimality - Hybrid of QR and DG -------------------
            [time_QD(CNT,w+1), H_QD, sensors_QD] = F_sensor_QD(U,p);
            det_QD (CNT,w+1) = F_calc_det  (p,H_QD,U);
            tr_QD  (CNT,w+1) = F_calc_trace(p,H_QD,U);
            eig_QD (CNT,w+1) = F_calc_eigen(p,H_QD,U);
            
            %% A-optimality -  Greedy -------------------------------
            [time_AG(CNT,w+1), H_AG, sensors_AG] = F_sensor_AG(U,p);
            det_AG (CNT,w+1) = F_calc_det  (p,H_AG,U);
            tr_AG  (CNT,w+1) = F_calc_trace(p,H_AG,U);
            eig_AG (CNT,w+1) = F_calc_eigen(p,H_AG,U);
            
            %% E-optimality -  Greedy -------------------------------
            [time_EG(CNT,w+1), H_EG, sensors_EG] = F_sensor_EG(U,p);
            det_EG (CNT,w+1) = F_calc_det  (p,H_EG,U);
            tr_EG  (CNT,w+1) = F_calc_trace(p,H_EG,U);
            eig_EG (CNT,w+1) = F_calc_eigen(p,H_EG,U);
        end
        
        %% Averaging ================================================
        [ time_rand, det_rand, tr_rand, eig_rand ]...
            = F_data_ave1( CNT, num_ave, time_rand, det_rand, tr_rand, eig_rand );
        [ time_DC, det_DC, tr_DC, eig_DC ]...
            = F_data_ave1( CNT, num_ave, time_DC, det_DC, tr_DC, eig_DC );
        [ time_QR, det_QR, tr_QR, eig_QR ]...
            = F_data_ave1( CNT, num_ave, time_QR, det_QR, tr_QR, eig_QR );
        [ time_DG, det_DG, tr_DG, eig_DG ]...
            = F_data_ave1( CNT, num_ave, time_DG, det_DG, tr_DG, eig_DG );
        [ time_QD, det_QD, tr_QD, eig_QD ]...
            = F_data_ave1( CNT, num_ave, time_QD, det_QD, tr_QD, eig_QD );
        [ time_AG, det_AG, tr_AG, eig_AG ]...
            = F_data_ave1( CNT, num_ave, time_AG, det_AG, tr_AG, eig_AG );
        [ time_EG, det_EG, tr_EG, eig_EG ]...
            = F_data_ave1( CNT, num_ave, time_EG, det_EG, tr_EG, eig_EG );
        NT_TOL_cal_DC(CNT,1)=mean(NT_TOL_cal_DC(CNT,2:w+1));
        iter_DC(CNT,1)=mean(iter_DC(CNT,2:w+1));
        
        %% Sensor location ==========================================
        sensor_memo = zeros(p,7);
        sensor_memo(1:p,1) = sensors_rand(1:p)';
        sensor_memo(1:p,2) = sensors_DC(1:p);
        sensor_memo(1:p,3) = sensors_QR(1:p)';
        sensor_memo(1:p,4) = sensors_DG(1:p);
        sensor_memo(1:p,5) = sensors_QD(1:p)';
        sensor_memo(1:p,6) = sensors_AG(1:p)';
        sensor_memo(1:p,7) = sensors_EG(1:p)';
        filename = [workdir, '/sensors_p_', num2str(p), '.mat'];
        save(filename,'sensor_memo');
        
        text = [ '---> ', num2str(p), ' sensor selection finished!' ];
        disp(text);
    end
end

%%NOAA-SST =========================================================
if num_problem == 2
    
    %% Preprocces for NOAA-SST ======================================
    text='Readinng/Arranging a NOAA-SST dataset';
    disp(text);
    [Lat, Lon, time, mask, sst]...
        = F_pre_read_NOAA_SST( ['sst.wkmean.1990-present.nc'], ['lsmask.nc'] );
    %  [Uorg, Sorg, Vorg, Xorg, meansst, n, Xtest] = F_pre_SVD_NOAA_SST_training_test(m, time, mask, sst,NumCrossVal,CV);
    
    [Uorg, Sorg, Vorg, Xorg, meansst, n] = F_pre_SVD_NOAA_SST(m, time, mask, sst);
    %    F_map_original(num_video, Xorg, meansst, mask, time, videodir);
    [U, Error_ave_pod, Error_std_pod]...
        = F_pre_truncatedSVD(r, Xorg, Uorg, Sorg, Vorg, num_video, meansst, mask, time, m, videodir);
    Error_ave_pod = repmat( Error_ave_pod , size(ps,2) );
    text='Complete Reading/Arranging a NOAA-SST dataset!';
    disp(text);
    
    det_CTC = det(U'*U);
    
    
    [Q,R,~] = qr(U','vector');
    for tt=1:n
        U_QR=Q*R(:,1:tt);
%         size(Q);
%         size(R);
        CTC=U_QR*U_QR';
        resi(tt,1)=norm(U'*U-CTC);
    end
    
    aaa
    %% QR
    sensor=[9660; 21899; 9426; 7000; 5563; 11316; 14807; 25102; 17970; 14230; 9663; 10793; 21607; 12459; 6701; 3764; 9953; 24813; 10151; 13997];
    for p=1:20
        [H]=F_calc_sensormatrix(p, n, sensor);
        det_QR (p,1) = F_calc_det_QR  (p,H,U);
    end
    
    %% DG
    sensor=[9660; 21899; 9426; 7000; 5563; 11316; 14807; 25102; 17970; 14230; 10792; 21278; 8848; 10349; 9664; 21903; 31372; 11089; 17393; 8387];
    for p=1:20
        [H]=F_calc_sensormatrix(p, n, sensor);
        det_DG (p,1) = F_calc_det_QR  (p,H,U);
    end
    
    
    
    aaa
    
    CNT=0;
    %% Sensor selection =============================================
    for p = ps
        CNT = CNT+1;
        text = [ num2str(p),' sensor selection started --->' ];
        disp(text);
        
        %% Random selection -----------------------------------------
        % Average loop
        for w=1:1:num_ave
            [time_rand(CNT,w+1), H_rand, sensors_rand] = F_sensor_random(n,p);
            det_rand(CNT,w+1) = F_calc_det  (p,H_rand,U);
            tr_rand (CNT,w+1) = F_calc_trace(p,H_rand,U);
            eig_rand(CNT,w+1) = F_calc_eigen(p,H_rand,U);
        end
        % Averaging
        %         [ time_rand, det_rand, tr_rand, eig_rand, Error_rand, Error_std_rand ]...
        %             = F_data_ave2( CNT, num_ave, time_rand, det_rand, tr_rand, eig_rand, Error_rand, Error_std_rand );
        
        %% D-optimality - Convex-------------------------------------
        %!! This is very time consuming proceduce, We do not recommend to try this
        %         [time_DC(CNT,1), H_DC, sensors_DC, ...
        %             NT_TOL_cal_DC(CNT,1), iter_DC(CNT,1)] ...
        %             = F_sensor_DC(U,p,maxiteration);
        %         det_DC(CNT,1) = F_calc_det(p,H_DC,U);
        %         tr_DC(CNT,1)  = F_calc_trace(p,H_DC,U);
        %         eig_DC(CNT,1) = F_calc_eigen(p,H_DC,U);
        %!! I recommend you use the following dummy values
        %   if you do not need the solution in the convex approximation in NOAA-SST.
        time_DC(CNT,1) = time_rand(CNT,1);
        det_DC (CNT,1) = det_rand (CNT,1);
        tr_DC  (CNT,1) = tr_rand  (CNT,1);
        eig_DC (CNT,1) = eig_rand (CNT,1);
        H_DC=H_rand;
        sensors_DC=sensors_rand;
        NT_TOL_cal_DC(CNT,w+1)=0;
        iter_DC(CNT,w+1)=0;
        %             NT_TOL_cal_DC(CNT,1)=mean(NT_TOL_cal_DC(CNT,2:w+1));
        %             iter_DC(CNT,1)=mean(iter_DC(CNT,2:w+1));
        
        %% Maximization of row norm - Greedy based on QR ------------
        [time_QR(CNT,1), H_QR, sensors_QR] = F_sensor_QR(U,p);
        det_QR (CNT,1) = F_calc_det_QR  (p,H_QR,U);
        tr_QR  (CNT,1) = F_calc_trace(p,H_QR,U);
        eig_QR (CNT,1) = F_calc_eigen(p,H_QR,U);
        
        %         time_QR(CNT,1) = time_rand(CNT,1);
        %         det_QR(CNT,1) = det_rand (CNT,1);
        %         tr_QR(CNT,1) = tr_rand  (CNT,1);
        %         eig_QR(CNT,1) = eig_rand (CNT,1);
        %         H_QR=H_rand;
        %         sensors_QR=sensors_rand;
        
        %% D-optimality - Greedy ------------------------------------
        [time_DG(CNT,1), H_DG, sensors_DG] = F_sensor_DG(U,p);
        det_DG (CNT,1) = F_calc_det  (p,H_DG,U);
        tr_DG  (CNT,1) = F_calc_trace(p,H_DG,U);
        eig_DG (CNT,1) = F_calc_eigen(p,H_DG,U);
        
        %% D-optimality - Hybrid of QR and DG -----------------------
        [time_QD(CNT,1), H_QD, sensors_QD] = F_sensor_QD(U,p);
        det_QD (CNT,1) = F_calc_det  (p,H_QD,U);
        tr_QD  (CNT,1) = F_calc_trace(p,H_QD,U);
        eig_QD (CNT,1) = F_calc_eigen(p,H_QD,U);
        
        %% A-optimality -  Greedy -----------------------------------
        [time_AG(CNT,1), H_AG, sensors_AG] = F_sensor_AG(U,p);
        det_AG (CNT,1) = F_calc_det  (p,H_AG,U);
        tr_AG  (CNT,1) = F_calc_trace(p,H_AG,U);
        eig_AG (CNT,1) = F_calc_eigen(p,H_AG,U);
        
        
        %% E-optimality -  Greedy -----------------------------------
        [time_EG(CNT,1), H_EG, sensors_EG] = F_sensor_EG(U,p);
        det_EG (CNT,1) = F_calc_det  (p,H_EG,U);
        tr_EG  (CNT,1) = F_calc_trace(p,H_EG,U);
        eig_EG (CNT,1) = F_calc_eigen(p,H_EG,U);
        
        
        %% GappyPOD+R ------------------
        for w=1:1:num_ave
            [time_GappyPODR(CNT,w+1), H_GappyPODR, sensors_GappyPODR] = F_sensor_GappyPODR(U, sensors_QR, p);
            det_GappyPODR(CNT,w+1) = F_calc_det  (p,H_GappyPODR,U);
            tr_GappyPODR(CNT,w+1) = F_calc_trace(p,H_GappyPODR,U);
            eig_GappyPODR(CNT,w+1) = F_calc_eigen(p,H_GappyPODR,U);
            
        end
        % Averaging
        %         [ time_GappyPODR, det_GappyPODR, tr_GappyPODR, eig_GappyPODR, Error_GappyPODR, Error_std_GappyPODR]...
        %             = F_data_ave2( CNT, num_ave,  time_GappyPODR, det_GappyPODR, tr_GappyPODR, eig_GappyPODR, Error_GappyPODR, Error_std_GappyPODR);
        %
        %% GappyPOD+E ------------------
        [time_GappyPODE(CNT,1), H_GappyPODE, sensors_GappyPODE] = F_sensor_GappyPODE(U,p);
        det_GappyPODE(CNT,1) = F_calc_det  (p,H_GappyPODE, U);
        tr_GappyPODE(CNT,1) = F_calc_trace(p,H_GappyPODE, U);
        eig_GappyPODE(CNT,1) = F_calc_eigen(p,H_GappyPODE,U);
        
        
        %% Sensor location ==========================================
        sensor_memo = zeros(p,9);
        sensor_memo(1:p,1) = sensors_rand(1:p)';
        sensor_memo(1:p,2) = sensors_DC(1:p);
        sensor_memo(1:p,3) = sensors_QR(1:p)';
        sensor_memo(1:p,4) = sensors_DG(1:p);
        sensor_memo(1:p,5) = sensors_QD(1:p)';
        sensor_memo(1:p,6) = sensors_AG(1:p)';
        sensor_memo(1:p,7) = sensors_EG(1:p)';
        sensor_memo(1:p,8) = sensors_GappyPODR(1:p)';
        sensor_memo(1:p,9) = sensors_GappyPODE(1:p)';
        
        filename = [workdir, '/sensors_p_', num2str(p), '.mat'];
        save(filename,'sensor_memo');
        
        %% Video ====================================================
        %             name='rand';
        %             F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
        %                 sensors_rand, Zestimate_rand, name, videodir, sensordir)
        %             name='DC';
        %             F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
        %                 sensors_DC, Zestimate_DC, name, videodir, sensordir)
        %             name='QR';
        %             F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
        %                 sensors_QR, Zestimate_QR, name, videodir, sensordir)
        %             name='DG';
        %             F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
        %                 sensors_DG, Zestimate_DG, name, videodir, sensordir)
        %             name='QD';
        %             F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
        %                 sensors_QD, Zestimate_QD, name, videodir, sensordir)
        %             name='AG';
        %             F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
        %                 sensors_AG, Zestimate_AG, name, videodir, sensordir)
        %             name='EG';
        %             F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
        %                 sensors_EG, Zestimate_EG, name, videodir, sensordir)
        %
        %             text = [ '---> ', num2str(p), ' sensor selection finished!' ];
        %             disp(text);
    end
    
    
    %% Data organization ================================================
    % Arrange
    [time_all] = F_data_arrange1( ps,   CNT, time_rand, time_DC, time_QR,...
        time_DG,   time_QD,   time_AG,     time_EG );
    [det_all]  = F_data_arrange1( ps,   CNT, det_rand,  det_DC,  det_QR, ...
        det_DG,    det_QD,    det_AG,      det_EG  );
    [tr_all]   = F_data_arrange1( ps,   CNT, tr_rand,   tr_DC,   tr_QR,...
        tr_DG,     tr_QD,     tr_AG,       tr_EG   );
    [eig_all]  = F_data_arrange1( ps,   CNT, eig_rand,  eig_DC,  eig_QR,...
        eig_DG,    eig_QD,    eig_AG,      eig_EG  );
    %     if num_problem == 2
    %         [Error] = F_data_arrange4( ps,         CNT, ...
    %             Error_rand(:,1), Error_std_rand(:,1), ...
    %             Error_DC,   Error_std_DC,   ...
    %             Error_QR,   Error_std_QR,   ...
    %             Error_DG,   Error_std_DG,   ...
    %             Error_QD,   Error_std_QD,   ...
    %             Error_AG,   Error_std_AG,   ...
    %             Error_EG,   Error_std_EG,   ...
    %             Error_GappyPODE,   Error_std_GappyPODE,   ...
    %             Error_GappyPODR(:,1),   Error_std_GappyPODR(:,1),   ...
    %             Error_ave_pod );
    %     end
    [log_DC] = F_data_arrange3( ps, CNT, NT_TOL_cal_DC, iter_DC );
    
    % Normalize
    [Normalized_det] = F_data_normalize( ps, CNT, det_rand, det_DC, det_QR, ...
        det_DG,  det_QD,   det_AG, det_EG );
    [Normalized_tr]  = F_data_normalize( ps, CNT, tr_rand,  tr_DC,  tr_QR,  ...
        tr_DG,   tr_QD,    tr_AG,  tr_EG  );
    [Normalized_eig] = F_data_normalize( ps, CNT, eig_rand, eig_DC, eig_QR, ...
        eig_DG,  eig_QD,   eig_AG, eig_EG );
    
    %% Save =============================================================
    cd(workdir)
    save('time.mat','time_all');
    save('det.mat','det_all');
    save('trace.mat','tr_all');
    save('eigen.mat','eig_all');
    save('Normalized_det.mat','Normalized_det');
    save('Normalized_trace.mat','Normalized_tr');
    save('Normalized_eigen.mat','Normalized_eig');
    save('time_rand.mat','time_rand');
    save('det_rand.mat','det_rand');
    save('trace_rand.mat','tr_rand');
    save('eigen_rand.mat','eig_rand');
    if num_problem == 1
        save('time_DC.mat','time_DC');
        save('time_QR.mat','time_QR');
        save('time_DG.mat','time_DG');
        save('time_QD.mat','time_QD');
        save('time_AG.mat','time_AG');
        save('time_EG.mat','time_EG');
        save('det_DC.mat','det_DC');
        save('det_QR.mat','det_QR');
        save('det_DG.mat','det_DG');
        save('det_QD.mat','det_QD');
        save('det_AG.mat','det_AG');
        save('det_EG.mat','det_EG');
    end
    save('log_DC.mat','log_DC');
    cd ../src
end



%% Cross-validation NOAA-SST =========================================================
if num_problem == 6 %    Cross-validation
    for NumCrossVal=1:CV
        time_NumCrossVal=tic;
        text = ['Cross Validation: ' num2str(NumCrossVal),' started *************************' ];
        disp(text);
        
        %% Preprocces for NOAA-SST ======================================
        text='Readinng/Arranging a NOAA-SST dataset';
        disp(text);
        [Lat, Lon, time, mask, sst]...
            = F_pre_read_NOAA_SST( ['sst.wkmean.1990-present.nc'], ['lsmask.nc'] );
        [Uorg, Sorg, Vorg, Xorg, meansst, n, Xtest] = F_pre_SVD_NOAA_SST_training_test(m, time, mask, sst,NumCrossVal,CV);
        
        U = Uorg(:,1:r);
        
        %         Xorg_test=U*pinv(U)*Xtest;
        %         [Utest,Stest,Vtest] = svd(Xtest,'econ');
        %         Xorg_test = Utest(:,1:r)*Stest(1:r,1:r)*Vtest(:,1:r)';
        %         Error(1:m/CV) = norm(Xorg_test(:,1:m/CV)-Xtest(:,1:m/CV)) / norm(Xtest(:,1:m/CV));
        %         Error_ave_pod = mean(Error);
        %         Error_std_pod = std(Error);
        %X = U*Sorg(1:r,1:r)*Vorg(:,1:r)';
        %      Error_ave_pod = repmat( Error_ave_pod , size(ps,2) );
        Error_ave_pod =0;
        Error_std_pod =0;
        text='Complete Reading/Arranging a NOAA-SST dataset!';
        disp(text);
        CNT=0;
        %% Sensor selection =============================================
        for p = ps
            CNT = CNT+1;
            text = [ num2str(p),' sensor selection started --->' ];
            disp(text);
            
            %% Random selection -----------------------------------------
            % Average loop
            for w=1:1:num_ave
                [time_rand(CNT,w+1), H_rand, sensors_rand] = F_sensor_random(n,p);
                det_rand(CNT,w+1) = F_calc_det  (p,H_rand,U);
                tr_rand (CNT,w+1) = F_calc_trace(p,H_rand,U);
                eig_rand(CNT,w+1) = F_calc_eigen(p,H_rand,U);
                [Zestimate_rand, Error_rand(CNT,w+1), Error_std_rand(CNT,w+1)] ...
                    = F_calc_error(m/CV, Xtest, U, H_rand); %Test
            end
            % Averaging
            [ time_rand, det_rand, tr_rand, eig_rand, Error_rand, Error_std_rand ]...
                = F_data_ave2( CNT, num_ave, time_rand, det_rand, tr_rand, eig_rand, Error_rand, Error_std_rand );
            
            %% D-optimality - Convex-------------------------------------
            %!! This is very time consuming proceduce, We do not recommend to try this
            [time_DC(CNT,1), H_DC, sensors_DC, ...
                NT_TOL_cal_DC(CNT,1), iter_DC(CNT,1)] ...
                = F_sensor_DC(U,p,maxiteration);
            det_DC(CNT,1) = F_calc_det(p,H_DC,U);
            tr_DC(CNT,1)  = F_calc_trace(p,H_DC,U);
            eig_DC(CNT,1) = F_calc_eigen(p,H_DC,U);
            %!! I recommend you use the following dummy values
            %   if you do not need the solution in the convex approximation in NOAA-SST.
            %             time_DC(CNT,1) = time_rand(CNT,1);
            %             det_DC (CNT,1) = det_rand (CNT,1);
            %             tr_DC  (CNT,1) = tr_rand  (CNT,1);
            %             eig_DC (CNT,1) = eig_rand (CNT,1);
            %             H_DC=H_rand;
            %             sensors_DC=sensors_rand;
            %             NT_TOL_cal_DC(CNT,w+1)=0;
            %             iter_DC(CNT,w+1)=0;
            %!!
            [Zestimate_DC, Error_DC(CNT,1), Error_std_DC(CNT,1)] ...
                = F_calc_error(m/CV, Xorg, U, H_DC);
            %             NT_TOL_cal_DC(CNT,1)=mean(NT_TOL_cal_DC(CNT,2:w+1));
            %             iter_DC(CNT,1)=mean(iter_DC(CNT,2:w+1));
            
            %% Maximization of row norm - Greedy based on QR ------------
            %             [time_QR(CNT,1), H_QR, sensors_QR] = F_sensor_QR(U,p);
            %             det_QR (CNT,1) = F_calc_det  (p,H_QR,U);
            %             tr_QR  (CNT,1) = F_calc_trace(p,H_QR,U);
            %             eig_QR (CNT,1) = F_calc_eigen(p,H_QR,U);
            
            time_QR(CNT,1) = time_rand(CNT,1);
            det_QR(CNT,1) = det_rand (CNT,1);
            tr_QR(CNT,1) = tr_rand  (CNT,1);
            eig_QR(CNT,1) = eig_rand (CNT,1);
            H_QR=H_rand;
            sensors_QR=sensors_rand;
            
            [Zestimate_QR, Error_QR(CNT,1), Error_std_QR(CNT,1)] ...
                = F_calc_error(m/CV, Xtest, U, H_QR);%Test
            
            %% D-optimality - Greedy ------------------------------------
            [time_DG(CNT,1), H_DG, sensors_DG] = F_sensor_DG(U,p);
            det_DG (CNT,1) = F_calc_det  (p,H_DG,U);
            tr_DG  (CNT,1) = F_calc_trace(p,H_DG,U);
            eig_DG (CNT,1) = F_calc_eigen(p,H_DG,U);
            [Zestimate_DG, Error_DG(CNT,1), Error_std_DG(CNT,1)] ...
                = F_calc_error(m/CV, Xtest, U, H_DG);%Test
            
            %% D-optimality - Hybrid of QR and DG -----------------------
            [time_QD(CNT,1), H_QD, sensors_QD] = F_sensor_QD(U,p);
            det_QD (CNT,1) = F_calc_det  (p,H_QD,U);
            tr_QD  (CNT,1) = F_calc_trace(p,H_QD,U);
            eig_QD (CNT,1) = F_calc_eigen(p,H_QD,U);
            [Zestimate_QD, Error_QD(CNT,1), Error_std_QD(CNT,1)] ...
                = F_calc_error(m/CV, Xtest, U, H_QD);%Test
            
            %% A-optimality -  Greedy -----------------------------------
            [time_AG(CNT,1), H_AG, sensors_AG] = F_sensor_AG(U,p);
            det_AG (CNT,1) = F_calc_det  (p,H_AG,U);
            tr_AG  (CNT,1) = F_calc_trace(p,H_AG,U);
            eig_AG (CNT,1) = F_calc_eigen(p,H_AG,U);
            [Zestimate_AG, Error_AG(CNT,1), Error_std_AG(CNT,1)] ...
                = F_calc_error(m/CV, Xtest, U, H_AG);%Test
            
            %% E-optimality -  Greedy -----------------------------------
            [time_EG(CNT,1), H_EG, sensors_EG] = F_sensor_EG(U,p);
            det_EG (CNT,1) = F_calc_det  (p,H_EG,U);
            tr_EG  (CNT,1) = F_calc_trace(p,H_EG,U);
            eig_EG (CNT,1) = F_calc_eigen(p,H_EG,U);
            [Zestimate_EG, Error_EG(CNT,1), Error_std_EG(CNT,1)] ...
                = F_calc_error(m/CV, Xtest, U, H_EG);%Test
            
            %% GappyPOD+R ------------------
            for w=1:1:num_ave
                [time_GappyPODR(CNT,w+1), H_GappyPODR, sensors_GappyPODR] = F_sensor_GappyPODR(U, sensors_QR, p);
                det_GappyPODR(CNT,w+1) = F_calc_det  (p,H_GappyPODR,U);
                tr_GappyPODR(CNT,w+1) = F_calc_trace(p,H_GappyPODR,U);
                eig_GappyPODR(CNT,w+1) = F_calc_eigen(p,H_GappyPODR,U);
                [Zestimate_GappyPODR, Error_GappyPODR(CNT,w+1), Error_std_GappyPODR(CNT,w+1)] ...
                    = F_calc_error(m/CV, Xtest, U, H_GappyPODR);%Test
            end
            % Averaging
            [ time_GappyPODR, det_GappyPODR, tr_GappyPODR, eig_GappyPODR, Error_GappyPODR, Error_std_GappyPODR]...
                = F_data_ave2( CNT, num_ave,  time_GappyPODR, det_GappyPODR, tr_GappyPODR, eig_GappyPODR, Error_GappyPODR, Error_std_GappyPODR);
            
            %% GappyPOD+E ------------------
            [time_GappyPODE(CNT,1), H_GappyPODE, sensors_GappyPODE] = F_sensor_GappyPODE(U,p);
            det_GappyPODE(CNT,1) = F_calc_det  (p,H_GappyPODE, U);
            tr_GappyPODE(CNT,1) = F_calc_trace(p,H_GappyPODE, U);
            eig_GappyPODE(CNT,1) = F_calc_eigen(p,H_GappyPODE,U);
            [Zestimate_GappyPODE, Error_GappyPODE(CNT,1), Error_std_GappyPODE(CNT,1)] ...
                = F_calc_error(m/CV, Xtest, U, H_GappyPODE);%Test
            
            %% Sensor location ==========================================
            sensor_memo = zeros(p,9);
            sensor_memo(1:p,1) = sensors_rand(1:p)';
            sensor_memo(1:p,2) = sensors_DC(1:p);
            sensor_memo(1:p,3) = sensors_QR(1:p)';
            sensor_memo(1:p,4) = sensors_DG(1:p);
            sensor_memo(1:p,5) = sensors_QD(1:p)';
            sensor_memo(1:p,6) = sensors_AG(1:p)';
            sensor_memo(1:p,7) = sensors_EG(1:p)';
            sensor_memo(1:p,8) = sensors_GappyPODR(1:p)';
            sensor_memo(1:p,9) = sensors_GappyPODE(1:p)';
            
            filename = [workdir, '/sensors_NumCV_', num2str(NumCrossVal), '_p_', num2str(p), '.mat'];
            save(filename,'sensor_memo');
            
            %% Video ====================================================
            %             name='rand';
            %             F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
            %                 sensors_rand, Zestimate_rand, name, videodir, sensordir)
            %             name='DC';
            %             F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
            %                 sensors_DC, Zestimate_DC, name, videodir, sensordir)
            %             name='QR';
            %             F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
            %                 sensors_QR, Zestimate_QR, name, videodir, sensordir)
            %             name='DG';
            %             F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
            %                 sensors_DG, Zestimate_DG, name, videodir, sensordir)
            %             name='QD';
            %             F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
            %                 sensors_QD, Zestimate_QD, name, videodir, sensordir)
            %             name='AG';
            %             F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
            %                 sensors_AG, Zestimate_AG, name, videodir, sensordir)
            %             name='EG';
            %             F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
            %                 sensors_EG, Zestimate_EG, name, videodir, sensordir)
            %
            %             text = [ '---> ', num2str(p), ' sensor selection finished!' ];
            %             disp(text);
        end
        
        
        %% Data organization ================================================
        % Arrange
        [time_all] = F_data_arrange1( ps,   CNT, time_rand, time_DC, time_QR,...
            time_DG,   time_QD,   time_AG,     time_EG );
        [det_all]  = F_data_arrange1( ps,   CNT, det_rand,  det_DC,  det_QR, ...
            det_DG,    det_QD,    det_AG,      det_EG  );
        [tr_all]   = F_data_arrange1( ps,   CNT, tr_rand,   tr_DC,   tr_QR,...
            tr_DG,     tr_QD,     tr_AG,       tr_EG   );
        [eig_all]  = F_data_arrange1( ps,   CNT, eig_rand,  eig_DC,  eig_QR,...
            eig_DG,    eig_QD,    eig_AG,      eig_EG  );
        if num_problem == 2
            [Error] = F_data_arrange4( ps,         CNT, ...
                Error_rand(:,1), Error_std_rand(:,1), ...
                Error_DC,   Error_std_DC,   ...
                Error_QR,   Error_std_QR,   ...
                Error_DG,   Error_std_DG,   ...
                Error_QD,   Error_std_QD,   ...
                Error_AG,   Error_std_AG,   ...
                Error_EG,   Error_std_EG,   ...
                Error_GappyPODE,   Error_std_GappyPODE,   ...
                Error_GappyPODR(:,1),   Error_std_GappyPODR(:,1),   ...
                Error_ave_pod );
        end
        [log_DC] = F_data_arrange3( ps, CNT, NT_TOL_cal_DC, iter_DC );
        
        % Normalize
        [Normalized_det] = F_data_normalize( ps, CNT, det_rand, det_DC, det_QR, ...
            det_DG,  det_QD,   det_AG, det_EG );
        [Normalized_tr]  = F_data_normalize( ps, CNT, tr_rand,  tr_DC,  tr_QR,  ...
            tr_DG,   tr_QD,    tr_AG,  tr_EG  );
        [Normalized_eig] = F_data_normalize( ps, CNT, eig_rand, eig_DC, eig_QR, ...
            eig_DG,  eig_QD,   eig_AG, eig_EG );
        
        %% Save =============================================================
        cd(workdir)
        save('time.mat','time_all');
        save('det.mat','det_all');
        save('trace.mat','tr_all');
        save('eigen.mat','eig_all');
        save('Normalized_det.mat','Normalized_det');
        save('Normalized_trace.mat','Normalized_tr');
        save('Normalized_eigen.mat','Normalized_eig');
        save('time_rand.mat','time_rand');
        save('det_rand.mat','det_rand');
        save('trace_rand.mat','tr_rand');
        save('eigen_rand.mat','eig_rand');
        if num_problem == 1
            save('time_DC.mat','time_DC');
            save('time_QR.mat','time_QR');
            save('time_DG.mat','time_DG');
            save('time_QD.mat','time_QD');
            save('time_AG.mat','time_AG');
            save('time_EG.mat','time_EG');
            save('det_DC.mat','det_DC');
            save('det_QR.mat','det_QR');
            save('det_DG.mat','det_DG');
            save('det_QD.mat','det_QD');
            save('det_AG.mat','det_AG');
            save('det_EG.mat','det_EG');
        end
        if num_problem == 2
            filename_error=['Error_CV_',num2str(NumCrossVal),'.mat'];
            save(filename_error,'Error');
            filename_error=['det_CV_',num2str(NumCrossVal),'.mat'];
            save(filename_error,'det_all');
            %   save('Error_rand.mat','Error_rand');
            
        end
        save('log_DC.mat','log_DC');
        cd ../src
        toc(time_NumCrossVal);
    end
end
warning('on','all')
disp('Congratulations!');
cd ../src
%% ///////////////////////////////////////////////////////////////////
%% Main program end