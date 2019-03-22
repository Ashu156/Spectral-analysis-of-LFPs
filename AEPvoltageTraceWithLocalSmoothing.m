close all;
clear all;
clc;

load('E:\LFPs\AEP2019\Rat27\22kSingleCall\matfile.mat');

% 1st time series (full)
ts_csc1 = (CSC26_TS)./10^6;
dp_csc1 = (CSC26_DP); % CSC data points in bit values
Fs = str2num(CSC26_NlxHeader{14}(end-3:end));  % Sampling frequency
Fnyq = round(Fs/2); % Nyquist frequency
ADBitVolts = str2num(CSC26_NlxHeader{15}(13:28)) ;% in Volts
ADV1 = 10^6*ADBitVolts; % in microVolts
dp_csc1 = dp_csc1*ADV1; % CSC data points in microVolts
R1 = dp_csc1(:); % Linearizing the data points
R1 = detrend(dp_csc1,'constant'); % detrending

tt = [0:1/Fs:511/Fs]' ;
tts_csc1 = [];
tts_csc2 = [];

for i = 1:length(ts_csc1)
    tts_csc1 = [tts_csc1  [ts_csc1(i) + tt]];
end

tts_csc1 = tts_csc1(:); % Linearized timestamps for the 1st time series

ts_events = (Events_tone_on + 000000.00)./10^6; % Selecting events of choice (onset of white noise pips in this case)

x_min = -0.5;
x_max = 1.0;
kmin = round(Fs*(x_min));         % Min x-limit  
kmax = round(Fs*(x_max));         % Max x-limit  

tx = (x_min:1/Fs:x_max); % time in milliseconds


numTrials = length(ts_events); % Number of trials

for i = 1:numTrials
    min_dist1(1,i) = min(abs(tts_csc1(:)-ts_events(i)));
    k1(1,i) = find(abs(tts_csc1(:)-ts_events(i)) == min_dist1(1,i)) ;
    eeg1(:,i) = R1(k1(1,i) + kmin : k1(1,i) + kmax);
end


F_cutL = 2; F_cutH = 12; %Low and High cut off frequencies (Hz)
[zz,pp,kk] = ellip(20, 0.2, 80, [F_cutL F_cutH]./Fnyq);
[sos,g] = zp2sos(zz,pp,kk);	      % Convert to SOS form
Hd = dfilt.df2tsos(sos,g);        % Create a dfilt object
filtered1  = filtfilthd(Hd, eeg1);  % Filtered signal
filtered1 = mean(filtered1,2);
%% change row to column
data = eeg1;
dtmp = [];
if isstruct(data);
   C = length(data);
   if C == 1;
      fnames = fieldnames(data);
      eval(['dtmp = data.' fnames{1} ';'])
      data = dtmp(:);
   end
else
  [N,C] = size(data);
  if N == 1 || C == 1;
    data = data(:);
  end;
end;

%% 

data = data';
win = [min(tx) max(tx)];
width = 10/Fs;
% err = 1;
T = win;
% T = [0 (N-1)/Fs];
% width = 50/Fs;
% err = 1;

t = min(T):1/Fs:max(T);

  indx = find(t>T(1) & t<T(2));
  t = t(indx);
  data = data(:,indx);


if width > (t(length(t))-t(1))/2
  disp('Width is too large for data segment: should be in seconds')
  disp('Turn off smoothing')
  width = 0;
end

s = t(2) - t(1); % sampling period
N = fix(width/s); % ratio of smoothing kernel to sampling period (how many sampling periods can be accommodated in the smoothing kernel)
NT = length(data(:,1)); % number of trials

if NT > 1 % if number of trials is more than one
    mdata = mean(data); % mean trace for all trials
else
    mdata = data; % mean trace is same as trial trace
end
if N > 4 % if number of sampling periods in a smoothing kernel is greater than 1, then, 
  smdata = locsmooth(mdata,N,fix(N/2)); % smoothened mean trace
% Fs = 1000; 
% Tw = size(mdata,2)/Fs; 
% Ts = Tw/2; 
% n = round(Fs*Tw);
% dn = round(Fs*Ts);
% mdata = mdata(:);
% nt = length(mdata);
% y_line = zeros(nt,1);
% norm = y_line;
% nwin = ceil((nt-n)/dn);
% yfit = zeros(nwin,n);
% xwt = ((1:n)- n/2)/(n/2);
% wt = (1-abs(xwt).^3).^3;
% for j = 1:nwin, 
% 	tseg = mdata(dn*(j-1)+1:dn*(j-1)+n);
% 	y1 = mean(tseg); 
% 	y2 = mean((1:n)'.*tseg)*2/(n+1);
% 	a = (y2-y1)*6/(n-1); 
%     b = y1-a*(n+1)/2;
% 	yfit(j,:) = (1:n)*a+b;
% 	y_line((j-1)*dn+(1:n)) = y_line((j-1)*dn+(1:n))+(yfit(j,:).*wt)';
% 	norm((j-1)*dn+(1:n)) = norm((j-1)*dn+(1:n))+wt';
% end

% mask = find(norm>0); 
% y_line(mask) = y_line(mask)./norm(mask);
% indx = (nwin-1)*dn+n-1;
% npts = length(mdata)-indx+1;
% y_line(indx:end) = (n+1:n+npts)'*a+b;
% tmp = y_line;
% % tmp2 = runline(data,n,dn); 
% smdata = tmp';
else
  smdata = mdata;  % smoothened trace is same as mean trace
end
  
% if errorbars requested then do a bootstrap over trials...

% Err = 0;
% if NT < 4; 
%   disp('Too few trials: no errorbars calculated')
%   err = 0;    
% end
% 
% if err ~= 0 && NT > 1
%   Nboot = 10;
%   bevk = 0;
%   sevk = 0;
%   for b = 1:Nboot
%     indx = floor(NT*rand(1,NT)) + 1;
%     evktmp = mean(data(indx,:));
%     if N > 4
%       evktmp = locsmooth(evktmp,N,fix(N/2));
%     end
%     bevk = bevk + evktmp;
%     sevk = sevk + evktmp.^2;
%   end
%   stdevk = sqrt((sevk/Nboot - bevk.^2/Nboot^2));
%   Err = stdevk;
% end

V = smdata; % smoothened mean voltage trace

plot(t,mdata,'k','linew',2) % plot mean voltage trace
hold on;
plot(t,smdata,'r','linew',2) % plot mean smoothened voltage trace
  
 mn1 = mean(mdata);  % mean of mean voltage trace
 mn2 = mean(smdata); % mean of smoothened mean voltage trace
 ax = get(gca,'xlim'); % get axis handle
 line(ax,mn1*[1 1],'Color','k') % plot a line equal to mean of mean voltage trace
 hold on
 line(ax,mn2*[0 0],'Color','r') % plot a line equal to mean of smoothened mean voltage trace
%   if err
%     line(ax,(mn+2*mean(stdevk))*[1 1],'color','r')
%     line(ax,(mn-2*mean(stdevk))*[1 1],'color','r')
%     hold off
%   end
legend('Raw voltage trace','Smoothened voltage trace');
%% local smoothening of voltage trace


