function [Psi,S,V,Xorg,Itest,mask, meansst, Y, x, time]=read_NOAA_SST(TT)
 text='Readinng/Arranging a NOAA-SST dataset';
 disp(text);

[ Lat, Lon, time, mask, sst ] = read_data_enso( ['sst.wkmean.1990-present.nc'],...
    ['lsmask.nc'] );
% each element of time array is a new week, in units of days
t0=datenum('1-Jan-1800 00:00:00')+time(1);
tfin=datenum('1-Jan-1800 00:00:00')+time(end);
datestr(t0,'yyyy/mm/dd');
datestr(tfin,'yyyy/mm/dd');

[m,n,p] = size(sst);
N = m*n;
X = zeros(N,p);

M = length(mask(mask==1));
Y = zeros(M,length(time));

for i=1:length(time)
    time_data(i)=datenum('1-Jan-1800 00:00:00')+time(i);
    snapshot = reshape(sst(:,:,i),N,1);
    Y(:,i) = snapshot(mask==1);
    % day_name(i)=datestr(time_data(i));
end


Iord = 1:length(time);
%Iord = 1:52*20;

%Itrain = Iord(1:52*16);
Itrain = Iord(1:TT); 
Itest=Itrain;
%Itest = Iord(~ismember(Iord,Itrain));

Train = Y(:,Itrain);
meansst = mean(Train,2);
Train = bsxfun(@minus,Train,meansst);

[Psi,S,V] = svd(Train,'econ');
Xorg=Psi*S*V';

for i=1:TT%num_video
    x(:,i) = Y(:,Itest(i))-meansst;
end

 text='Complete Reading/Arranging a NOAA-SST dataset!';
 disp(text);
end 