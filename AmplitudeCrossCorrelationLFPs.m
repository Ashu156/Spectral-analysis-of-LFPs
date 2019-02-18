% This piece of code calculates the cross-correlation between LFPs of two channels
% and plots the same. This is inspired by the following publication:
% "Cross-correlation of instantaneous amplitudes of field potential 
% oscillations: A straightforward method to estimate the directionality and lag between brain areas"
% https://doi.org/10.1016/j.jneumeth.2010.06.019

close all; 
clear all; 
clc;

load('F:\LFPs\AEP2019\Rat32\10kpips\matfile.mat');

% 1st time series (full)
ts_csc1 = (CSC28_TS)./10^6;
dp_csc1 = (CSC28_DP); %CSC data points in bit values
Fs = str2num(CSC28_NlxHeader{14}(end-3:end)); Fnyq = round(Fs/2); %sampling and nyquist frequency
ADBitVolts = str2num(CSC28_NlxHeader{15}(13:28)) ;% in Volts
ADV1 = 10^6*ADBitVolts; % in microVolts
dp_csc1 = dp_csc1*ADV1; %CSC data points in microVolts
R1 = dp_csc1(:); % Linearizing the data points
R1 = detrend(dp_csc1,'constant'); % linear detrending the 1st time series

% 2nd time series (full)
ts_csc2 = (CSC31_TS)./10^6;
dp_csc2 = (CSC31_DP); %CSC data points in bit values
Fs = str2num(CSC31_NlxHeader{14}(end-3:end)) ; Fnyq = round(Fs/2); %sampling and nyquist frequency
ADBitVolts = str2num(CSC31_NlxHeader{15}(13:38)) ;% in Volts
ADV2 = 10^6*ADBitVolts; % in microVolts
dp_csc1 = dp_csc2*ADV2; %CSC data points in microVolts
R2 = dp_csc2(:); % Linearizing the data points
R2 = detrend(dp_csc2,'constant'); % linear detrending the 2nd time series

% Linearized timestamps for the 1st and 2nd time series

tt = [0:1/Fs:511/Fs]' ;
tts_csc1 = [];
tts_csc2 = [];

for i = 1:length(ts_csc1)
    tts_csc1 = [tts_csc1  [ts_csc1(i) + tt]];
    tts_csc2 = [tts_csc2  [ts_csc2(i) + tt]];
end

tts_csc1 = tts_csc1(:); % Linearized timestamps for the 1st time series
tts_csc2 = tts_csc2(:); % Linearized timestamps for the 2nd time series


ts_events = (Events_tone_on + 000000.00)./10^6; % Selecting events of choice (onset of white noise pips in this case)

x_min = -0.2;
x_max = 1.0;
kmin = round(Fs*(x_min));         % Min x-limit  
kmax = round(Fs*(x_max));         % Max x-limit  

tx = (x_min:1/Fs:x_max)*1000 ; % time in milliseconds


numTrials = length(ts_events); % Number of trials

for i = 1:numTrials
    min_dist1(1,i) = min(abs(tts_csc1(:)-ts_events(i)));
    min_dist2(1,i) = min(abs(tts_csc2(:)-ts_events(i)));
    k1(1,i) = find(abs(tts_csc1(:)-ts_events(i)) == min_dist1(1,i)) ;
    k2(1,i) = find(abs(tts_csc2(:)-ts_events(i)) == min_dist2(1,i)) ;
    eeg1(:,i) = R1(k1(1,i) + kmin : k1(1,i) + kmax);
    eeg2(:,i) = R2(k2(1,i) + kmin : k2(1,i) + kmax);
end


F_cutL = 4; F_cutH = 10; %Low and High cut off frequencies (Hz)
[zz,pp,kk] = ellip(20, 0.2, 80, [F_cutL F_cutH]./Fnyq);
[sos,g] = zp2sos(zz,pp,kk);	      % Convert to SOS form
Hd = dfilt.df2tsos(sos,g);        % Create a dfilt object
filtered1  = filtfilthd(Hd, eeg1);  % Filtered signal
filtered1 = mean(filtered1,2);
filtered2  = filtfilthd(Hd, eeg2);  % Filtered signal
filtered2 = mean(filtered2,2);

filt_hilb1 = hilbert(filtered1); % calculating the Hilbert transform of eeg1
amp1 = abs(filt_hilb1); % calculating instantaneous amplitude of filtered eeg1

filt_hilb2 = hilbert(filtered2);    % calculating the Hilbert transform of eeg2
amp2 = abs(filt_hilb2);       % calculating instantaneous amplitude of filtered eeg2

% for i = 1:numTrials
%     amp12(:,i) = amp1(:,i) - mean(amp1(:,i));
%     amp22(:,i) = amp2(:,i) - mean(amp2(:,i));
% end

% for i =1:numTrials
%     [crosscorr(:,i), lags] = xcorr(amp1(:,i), amp2(:,i), 500, 'coeff'); % calculates crosscorrelation between amplitude vectors
%     lags(i,:) = (lags./Fs)*1e3; % converts lags into ms
%     g(i) = find(crosscorr(:,i) == max(crosscorr(:,i))); % identifies the index where the crosscorrelation peaks
%     max_crosscorr_lag(i) = lags(g(i)); % identifies the lag at which the crosscorrelation peaks
% end

[crosscorr, lags] = xcorr(amp1,amp2, 100, 'coeff');
lags = (lags./Fs)*1e3;
g = find(crosscorr == max(crosscorr));
max_crosscorr_lag = lags(g);

% lags = lags(1,:);

% for i = 1:numTrials
%     plot(lags,crosscorr(:,i));
% hold on
% end


figure('Color',[1,1,1]);
plot(lags,crosscorr','Color',[0 0 1],'Linewidth',2); hold on; % plots crosscorrelations as a function of lags
plot(lags(g),crosscorr(g),'rp','MarkerFaceColor',[1 0 0],'MarkerSize',10); % plots marker at the peak of the crosscorrelation
plot([0 0],[1.05*max(crosscorr) 0.95*min(crosscorr)],'Color',[0 0 0],'Linestyle',':','Linewidth',2); % plots dashed line at zero lag
set(gca,'XTick',[-100 -50 0 50 100]);
% axis tight; box off; xlim([-101 100]);
% xlabel('Lag(ms)','FontSize',14);
% ylabel('Crosscorrelation','FontSize',14);



%% end of code
