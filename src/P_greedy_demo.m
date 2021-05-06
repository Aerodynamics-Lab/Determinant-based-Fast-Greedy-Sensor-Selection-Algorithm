%% Main program
%% ///////////////////////////////////////////////////////////////////
% Comments:
% 	Collaborator: Yuji Saito, Keigo Yamada, Taku Nonomura
%                 Kumi Nakai, Takayuki Nagata
% 	Last modified: 2021/4/28
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
%num_problem=1; % //Randomized sensor problem//
%num_problem=2; % //NOAA-SST//
num_problem=3;% // Cross-validation NOAA-SST
% !<NOAA-SST> It takes a long time to obtain the solution in the convex
% !<NOAA-SST> approximation method and the convex method is commented out
% !<NOAA_SST> as default setting for reduction of demo time.
%
%% Parameters =======================================================
r = 10;
pmin = 1;
pinc = 1;
pmax = 20;
ps   = pmin:pinc:pmax;
CNT = 0; % Counter
maxiteration = 200; % Max iteration for convex approximation

%/////////////////////////////
% //Randomized sensor problem//
n = 2000;
num_ave = 1; % Number of iteration for averaging operation
%/////////////////////////////
% //NOAA-SST//
m = 52*10; % 10-years (52weeks/year)
num_video = 1; % maxmum: m
%/////////////////////////////
% // Cross-validation NOAA-SST//
CV=5;
%/////////////////////////////
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
            
            %% D-optimality - Convex---------------------------------
            %!! This is very time consuming proceduce, We do not recommend to try this
            [time_DC(CNT,w+1), H_DC, sensors_DC, ...
                NT_TOL_cal_DC(CNT,w+1), iter_DC(CNT,w+1)] ...
                = F_sensor_DC(U,p,maxiteration);
            det_DC (CNT,w+1) = F_calc_det  (p,H_DC,U);
            %!! I recommend you use the following dummy values
            %   if you do not need the solution in the convex approximation in NOAA-SST.
            %             time_DC(CNT,w+1) = time_rand(CNT,w+1);
            %             det_DC (CNT,w+1) = det_rand (CNT,w+1);
            %            H_DC=H_rand;
            %            sensors_DC=sensors_rand;
            %            NT_TOL_cal_DC(CNT,w+1)=0;
            %            iter_DC(CNT,w+1)=0;
            
            %% Maximization of row norm - Greedy based on QR --------
            [time_QR(CNT,w+1), H_QR, sensors_QR] = F_sensor_QR(U,p);
            det_QR (CNT,w+1) = F_calc_det  (p,H_QR,U);
            
            %% D-optimality - Greedy --------------------------------
            [time_DG(CNT,w+1), H_DG, sensors_DG] = F_sensor_DG(U,p);
            det_DG (CNT,w+1) = F_calc_det  (p,H_DG,U);
            
            %% D-optimality - Hybrid of QR and DG -------------------
            [time_QD(CNT,w+1), H_QD, sensors_QD] = F_sensor_QD(U,p);
            det_QD (CNT,w+1) = F_calc_det  (p,H_QD,U);
            
        end
        
        %% Averaging ================================================
        [ time_rand, det_rand]...
            = F_data_ave1( CNT, num_ave, time_rand, det_rand);
        [ time_DC, det_DC]...
            = F_data_ave1( CNT, num_ave, time_DC, det_DC);
        [ time_QR, det_QR]...
            = F_data_ave1( CNT, num_ave, time_QR, det_QR);
        [ time_DG, det_DG]...
            = F_data_ave1( CNT, num_ave, time_DG, det_DG);
        [ time_QD, det_QD]...
            = F_data_ave1( CNT, num_ave, time_QD, det_QD);
        NT_TOL_cal_DC(CNT,1)=mean(NT_TOL_cal_DC(CNT,2:w+1));
        iter_DC(CNT,1)=mean(iter_DC(CNT,2:w+1));
        
        %% Sensor location ==========================================
        sensor_memo = zeros(p,5);
        sensor_memo(1:p,1) = sensors_rand(1:p)';
        sensor_memo(1:p,2) = sensors_DC(1:p);
        sensor_memo(1:p,3) = sensors_QR(1:p)';
        sensor_memo(1:p,4) = sensors_DG(1:p);
        sensor_memo(1:p,5) = sensors_QD(1:p)';
        filename = [workdir, '/sensors_p_', num2str(p), '.mat'];
        save(filename,'sensor_memo');
        
        text = [ '---> ', num2str(p), ' sensor selection finished!' ];
        disp(text);
    end
    [time_all] = F_data_arrange1( ps,   CNT, time_rand, time_DC, time_QR,...
        time_DG,   time_QD);
    [det_all]  = F_data_arrange1( ps,   CNT, det_rand,  det_DC,  det_QR, ...
        det_DG,    det_QD);
    % Normalize
    [Normalized_det] = F_data_normalize( ps, CNT, det_rand, det_DC, det_QR, ...
        det_DG,  det_QD);
    
    cd(workdir)
    save('time.mat','time_all');
    save('det.mat','det_all');
    save('Normalized_det.mat','Normalized_det');
    save('time_rand.mat','time_rand');
    save('det_rand.mat','det_rand');
    save('time_DC.mat','time_DC');
    save('time_QR.mat','time_QR');
    save('time_DG.mat','time_DG');
    save('time_QD.mat','time_QD');
    save('det_DC.mat','det_DC');
    save('det_QR.mat','det_QR');
    save('det_DG.mat','det_DG');
    save('det_QD.mat','det_QD');
    
end

%%NOAA-SST =========================================================
if num_problem == 2
    
    %% Preprocces for NOAA-SST ======================================
    text='Readinng/Arranging a NOAA-SST dataset';
    disp(text);
    [Lat, Lon, time, mask, sst]...
        = F_pre_read_NOAA_SST( 'sst.wkmean.1990-present.nc', 'lsmask.nc');
    [Uorg, Sorg, Vorg, Xorg, meansst, n] = F_pre_SVD_NOAA_SST(m, time, mask, sst);
    [U, Error_ave_pod, Error_std_pod]...
        = F_pre_truncatedSVD(r, Xorg, Uorg, Sorg, Vorg, num_video, meansst, mask, time, m, videodir);
    %Error_ave_pod = repmat( Error_ave_pod , size(ps,2) );
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
            [Zestimate_rand, Error_rand(CNT,w+1), Error_std_rand(CNT,w+1)] ...
                = F_calc_error(m, Xorg, U, H_rand);
        end
        % Averaging
        [ time_rand, det_rand, Error_rand, Error_std_rand ]...
            = F_data_ave2(CNT, num_ave, time_rand, det_rand, Error_rand, Error_std_rand );
        
        %% D-optimality - Convex-------------------------------------
        %!! This is very time consuming proceduce, We do not recommend to try this
        %         [time_DC(CNT,1), H_DC, sensors_DC, ...
        %             NT_TOL_cal_DC(CNT,1), iter_DC(CNT,1)] ...
        %             = F_sensor_DC(U,p,maxiteration);
        %         det_DC(CNT,1) = F_calc_det(p,H_DC,U);
        
        %!! I recommend you use the following dummy values
        %   if you do not need the solution in the convex approximation in NOAA-SST.
        time_DC(CNT,1) = time_rand(CNT,1);
        det_DC (CNT,1) = det_rand (CNT,1);
        H_DC=H_rand;
        sensors_DC=sensors_rand;
        NT_TOL_cal_DC(CNT,w+1)=0;
        iter_DC(CNT,w+1)=0;
        
        [Zestimate_DC, Error_DC(CNT,1), Error_std_DC(CNT,1)] ...
            = F_calc_error(m, Xorg, U, H_DC);
        
        %% Maximization of row norm - Greedy based on QR ------------
        %                     [time_QR(CNT,1), H_QR, sensors_QR] = F_sensor_QR(U,p);
        %                     det_QR (CNT,1) = F_calc_det_QR  (p,H_QR,U);
        %!! I recommend you use the following dummy values in the oversampling cases
        time_QR(CNT,1) = time_rand(CNT,1);
        det_QR(CNT,1) = det_rand (CNT,1);
        H_QR=H_rand;
        sensors_QR=sensors_rand;
        
        [Zestimate_QR, Error_QR(CNT,1), Error_std_QR(CNT,1)] ...
            = F_calc_error(m, Xorg, U, H_QR);
        
        %% D-optimality - Greedy ------------------------------------
        [time_DG(CNT,1), H_DG, sensors_DG] = F_sensor_DG(U,p);
        det_DG (CNT,1) = F_calc_det  (p,H_DG,U);
        [Zestimate_DG, Error_DG(CNT,1), Error_std_DG(CNT,1)] ...
            = F_calc_error(m, Xorg, U, H_DG);
        
        %% D-optimality - Hybrid of QR and DG -----------------------
        [time_QD(CNT,1), H_QD, sensors_QD] = F_sensor_QD(U,p);
        det_QD (CNT,1) = F_calc_det  (p,H_QD,U);
        [Zestimate_QD, Error_QD(CNT,1), Error_std_QD(CNT,1)] ...
            = F_calc_error(m, Xorg, U, H_QD);
             
        %% Sensor location ==========================================
        sensor_memo = zeros(p,5);
        sensor_memo(1:p,1) = sensors_rand(1:p)';
        sensor_memo(1:p,2) = sensors_DC(1:p);
        sensor_memo(1:p,3) = sensors_QR(1:p)';
        sensor_memo(1:p,4) = sensors_DG(1:p);
        sensor_memo(1:p,5) = sensors_QD(1:p)';
        
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
        %
        %             text = [ '---> ', num2str(p), ' sensor selection finished!' ];
        %             disp(text);
    end
        
    %% Data organization ================================================
    % Arrange
    [time_all] = F_data_arrange1( ps,   CNT, time_rand, time_DC, time_QR,...
        time_DG,   time_QD);
    [det_all]  = F_data_arrange1( ps,   CNT, det_rand,  det_DC,  det_QR, ...
        det_DG,    det_QD);
    
    [Error] = F_data_arrange4( ps,         CNT, ...
        Error_rand(:,1), Error_std_rand(:,1), ...
        Error_DC,   Error_std_DC,   ...
        Error_QR,   Error_std_QR,   ...
        Error_DG,   Error_std_DG,   ...
        Error_QD,   Error_std_QD,   ...
        Error_ave_pod );
    [log_DC] = F_data_arrange3( ps, CNT, NT_TOL_cal_DC, iter_DC );
    
    
    %% Save =============================================================
    cd(workdir)
    save('time.mat','time_all');
    save('det.mat','det_all');
    save('error.mat','Error');
    save('time_rand.mat','time_rand');
    save('det_rand.mat','det_rand');
    save('log_DC.mat','log_DC');
    cd ../src
end

%% Cross-validation NOAA-SST =========================================================
if num_problem == 3 %    Cross-validation
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
                [time_rand(CNT,w+CV), H_rand, sensors_rand] = F_sensor_random(n,p);
                det_rand(CNT,w+CV) = F_calc_det  (p,H_rand,U);
                [Zestimate_rand, Error_rand(CNT,w+CV), Error_std_rand(CNT,w+CV)] ...
                    = F_calc_error(m/CV, Xtest, U, H_rand); %Test
            end
            % Averaging
            [time_rand, det_rand, Error_rand, Error_std_rand ]...
                = F_data_ave4(CNT, num_ave, CV, NumCrossVal, time_rand, det_rand,  Error_rand, Error_std_rand );
            
            %% D-optimality - Convex-------------------------------------
            %!! This is very time consuming proceduce, We do not recommend to try this
%             [time_DC(CNT,NumCrossVal+1), H_DC, sensors_DC, ...
%                 NT_TOL_cal_DC(CNT,1), iter_DC(CNT,1)] ...
%                 = F_sensor_DC(U,p,maxiteration);
%             det_DC(CNT,NumCrossVal+1) = F_calc_det(p,H_DC,U);
            %!! I recommend you use the following dummy values
            %   if you do not need the solution in the convex approximation in NOAA-SST.
                        time_DC(CNT,NumCrossVal+1) = time_rand(CNT,NumCrossVal+1);
                        det_DC (CNT,NumCrossVal+1) = det_rand (CNT,NumCrossVal+1);
                        H_DC=H_rand;
                        sensors_DC=sensors_rand;
                        NT_TOL_cal_DC(CNT,NumCrossVal+1)=0;
                        iter_DC(CNT,NumCrossVal+1)=0;
            %!!
            [Zestimate_DC, Error_DC(CNT,NumCrossVal+1), Error_std_DC(CNT,NumCrossVal+1)] ...
                = F_calc_error(m/CV, Xorg, U, H_DC);
            %             NT_TOL_cal_DC(CNT,1)=mean(NT_TOL_cal_DC(CNT,2:w+1));
            %             iter_DC(CNT,1)=mean(iter_DC(CNT,2:w+1));
            
            %% Maximization of row norm - Greedy based on QR ------------
            %             [time_QR(CNT,1), H_QR, sensors_QR] = F_sensor_QR(U,p);
            %             det_QR (CNT,1) = F_calc_det  (p,H_QR,U);
            %!! I recommend you use the following dummy values in the oversampling cases
            time_QR(CNT,NumCrossVal+1) = time_rand(CNT,1);
            det_QR(CNT,NumCrossVal+1) = det_rand (CNT,1);
            H_QR=H_rand;
            sensors_QR=sensors_rand;
            
            [Zestimate_QR, Error_QR(CNT,NumCrossVal+1), Error_std_QR(CNT,NumCrossVal+1)] ...
                = F_calc_error(m/CV, Xtest, U, H_QR);%Test
            
            %% D-optimality - Greedy ------------------------------------
            [time_DG(CNT,NumCrossVal+1), H_DG, sensors_DG] = F_sensor_DG(U,p);
            det_DG(CNT,NumCrossVal+1) = F_calc_det  (p,H_DG,U);
            [Zestimate_DG, Error_DG(CNT,NumCrossVal+1), Error_std_DG(CNT,NumCrossVal+1)] ...
                = F_calc_error(m/CV, Xtest, U, H_DG);%Test
            
            %% D-optimality - Hybrid of QR and DG -----------------------
            [time_QD(CNT,NumCrossVal+1), H_QD, sensors_QD] = F_sensor_QD(U,p);
            det_QD(CNT,NumCrossVal+1) = F_calc_det  (p,H_QD,U);
            [Zestimate_QD, Error_QD(CNT,NumCrossVal+1), Error_std_QD(CNT,NumCrossVal+1)] ...
                = F_calc_error(m/CV, Xtest, U, H_QD);%Test
            
            %% Sensor location ==========================================
            sensor_memo = zeros(p,5);
            sensor_memo(1:p,1) = sensors_rand(1:p)';
            sensor_memo(1:p,2) = sensors_DC(1:p);
            sensor_memo(1:p,3) = sensors_QR(1:p)';
            sensor_memo(1:p,4) = sensors_DG(1:p);
            sensor_memo(1:p,5) = sensors_QD(1:p)';
            
            filename = [workdir, '/sensors_NumCV_', num2str(NumCrossVal), '_p_', num2str(p), '.mat'];
            save(filename,'sensor_memo');
        end
        toc(time_NumCrossVal);
    end
    % averaging
     CNT=0;
     for p = ps
            CNT = CNT+1;
    [ time_rand, det_rand, Error_rand, Error_std_rand]...
        = F_data_ave3( CNT, NumCrossVal, time_rand, det_rand, Error_rand, Error_std_rand);
    [ time_DC, det_DC, Error_DC, Error_std_DC]...
        = F_data_ave3( CNT, NumCrossVal, time_DC, det_DC, Error_DC, Error_std_DC);
    [ time_QR, det_QR, Error_QR, Error_std_QR]...
        = F_data_ave3( CNT, NumCrossVal, time_QR, det_QR, Error_QR, Error_std_QR);
    [ time_DG, det_DG, Error_DG, Error_std_DG]...
        = F_data_ave3( CNT, NumCrossVal, time_DG, det_DG, Error_DG, Error_std_DG);
    [ time_QD, det_QD, Error_QD, Error_std_QD]...
        = F_data_ave3( CNT, NumCrossVal, time_QD, det_QD, Error_QD, Error_std_QD);
     end
    %% Data organization ================================================
    % Arranging
    [time_all] = F_data_arrange1( ps,   CNT, time_rand, time_DC, time_QR,...
        time_DG,   time_QD);
    [det_all]  = F_data_arrange1( ps,   CNT, det_rand,  det_DC,  det_QR, ...
        det_DG,    det_QD);
    
    [Error] = F_data_arrange5( ps,         CNT, ...
        Error_rand, Error_std_rand, ...
        Error_DC,   Error_std_DC,   ...
        Error_QR,   Error_std_QR,   ...
        Error_DG,   Error_std_DG,   ...
        Error_QD,   Error_std_QD,   ...
        Error_ave_pod );
    [log_DC] = F_data_arrange3( ps, CNT, NT_TOL_cal_DC, iter_DC );
    
    %% Save =============================================================
    cd(workdir)
    filename_error=['Time_CV_',num2str(NumCrossVal),'.mat'];
    save(filename_error,'time_all');
    filename_error=['Error_CV_',num2str(NumCrossVal),'.mat'];
    save(filename_error,'Error');
    filename_error=['det_CV_',num2str(NumCrossVal),'.mat'];
    save(filename_error,'det_all');
    save('log_DC.mat','log_DC');
    cd ../src
end
warning('on','all')
disp('Congratulations!');
cd ../src
%% ///////////////////////////////////////////////////////////////////
%% Main program end