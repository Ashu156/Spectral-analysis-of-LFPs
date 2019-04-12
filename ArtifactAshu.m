% This script uses discrete stationary wavelet decomposition of a signal
% with artifacts and outputs the cleaned signal. It is basically based on 
% wavelet denoising. The algorithm is described on
% the following publication:
% https://www.sciencedirect.com/science/article/pii/S0165027014000387?via%3Dihub
% https://www.ncbi.nlm.nih.gov/pubmed/24512692
% data_segment/data_art is the input noisy signal
% data_new is the cleaned output signal
% Written and tested in MATLAB 2014a
%% 
tic;       % Start a timer

clear all; % Clear the MATLAB workspace
close all; % Close any open window
clc;       % Clear the command window

%% Example dataset from paper

% Uncomment the below section to run the code on the example dataset
% provided in the paper. Make sure you have the dataset in your current
% workinfg directory. If not, then, make sure that you import the dataset 
% correctly.

% load('C:\Users\Dell\Downloads\ArtifactRejection\Array2_1_ch_1.mat');
% Fs = 40e3;
% data_length = 1e6; start_index = 1e3;
% data_segment = 10*ch(start_index + 1 : start_index + data_length)/2^15;
% 
% t = (0:length(data_segment)-1)./Fs;
% figure;hold on;
% plot(t, data_segment);
% xlabel('Time, Sec'); ylabel('Amplitude, mV');
% title('A Segment of Raw Data from Neural Recording');

%% Dataset from own recordings

[fileName, path] = uigetfile;
rawData = load(strcat(path, fileName));


x_min = -0.5; x_max = 1.5;

% CHOOSING THE REQUIRED CHANNEL

ts_csc = rawData.CSC26_TS; % Channel of choice
dp_csc = rawData.CSC26_DP; % Data points in channel of choice
Fs = str2double(rawData.CSC01_NlxHeader{14}(end-4:end)); % Sampling frequency
ADBitVolts = (str2double(rawData.CSC26_NlxHeader{15}(13:38)))*10^6;% in microVolts
dp_csc = dp_csc*ADBitVolts; %CSC data points in microVolts
dp_csc = dp_csc(:);

dp_csc = detrend(dp_csc,'constant');
ts_events = (rawData.Events_tone_on + 000000.00)./10^6; % Selecting events of choice (onset of white noise pips in this case)
numTrials = length(rawData.ts_events); % number of events
ts_csc = ts_csc./10^6 ; % In seconds
data_length = length(dp_csc);
data_segment = dp_csc(1:data_length);


% Generate all the timestamps corresponding to all data points
tt = [0:1/Fs:511/Fs]';
tts_csc = [];
for i = 1:length(ts_csc)
    tts_csc = [tts_csc  [ts_csc(i) + tt]];
end
tts_csc = tts_csc(:);

%% Initial Filtering and Threshold Calculation

data_art = data_segment;

data_length = 2^nextpow2 (length(data_art)); % length of the input signal

% zero-padding data for fast computation

if (data_length > length(data_art))
    data_length = 2^(nextpow2 (length(data_art))-1);
end

sig_us = data_art(1:data_length) - mean(data_art(1:data_length)); % signal detrending

Fnyq = round(Fs/2);
F_cutL = 150; F_cutH = 300; %Low and High cut off frequencies (Hz)
[zz,pp,kk] = ellip(20, 0.2, 80, [F_cutL F_cutH]./Fnyq);
[sos,g] = zp2sos(zz,pp,kk);	      % Convert to SOS form
Hd = dfilt.df2tsos(sos,g);        % Create a dfilt object
Data_filt  = filtfilthd(Hd, sig_us); % Filtered signal
x_bpf = Data_filt;


% x_bpf = bandpass_filter(sig_us, Fs, 150, 300, 512); % bandpass filtering the signal between 150 and 300 Hz



% x_bpf_sp = bandpass_filter(sig_us, Fs, 300, 6e3, 512);  % bandpass filtering the signal between 300 and 6000 Hz

reset(Hd);
F_cutL = 300; F_cutH = 490; %Low and High cut off frequencies (Hz)
[zz,pp,kk] = ellip(20, 0.2, 80, [F_cutL F_cutH]./Fnyq);
[sos,g] = zp2sos(zz,pp,kk);	      % Convert to SOS form
Hd = dfilt.df2tsos(sos,g);        % Create a dfilt object
Data_filt  = filtfilthd(Hd, sig_us); % Filtered signal
x_bpf_sp = Data_filt;

sp = median(abs(x_bpf_sp))/0.6745; % unbiased estimate of noise variance for spike data

thrp = sp*sqrt(2*log10(length(x_bpf_sp))); % determining universal threshold for spike data

s = median(abs(x_bpf))/0.6745; % estimating noise variance for LFP data

thr = s*sqrt(2*log10(length(x_bpf))); % determining universal threshold for LFP data

% x_hpf = highpass_filter(sig_us, Fs, 5e3, 512); % highpass filtering the signal above 5000 Hz

% sh = median(abs(x_hpf))/0.6745; % unbiased estimate of noise variance for highpass signal

% thrh = sh*sqrt(2*log10(length(x_hpf))); % determining universal threshold for highpass signal

%% Discrete Stationary Wavelet Transform(DSWT)

tic; % initialize timer

N = 10; % number of decompositional levels to be used

wave_name = 'haar'; % wavelet family to be used

%%--DSWT
[A,D] = swt(sig_us,N,wave_name); % performing DSWT

[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters(wave_name); % wavelet filter

%% Approximate Coefficient for Thresholding the signal

A_new = A(end,:); A_old = A(end,:); clear A;

min_ratio = min((A_new))/(median(abs(A_new))/0.6745);

max_ratio = max((A_new))/(median(abs(A_new))/0.6745);

avg_ratio = max(abs(A_new));

thr_ratio = 3;

if ( avg_ratio > 2*thr_ratio )
    k1 = 0.5;
    elseif ( 2*thr_ratio > avg_ratio > thr_ratio )
        k1 = 0.75;
else
    k1 = 1;
end

sigma = median(abs(A_new))/0.6745;
T = k1*sqrt(2*log10(length(A_new))*sigma^2); % modified threshold value
id = find((abs(A_new)> T)==1);
tau2 = 0.8;
lamda2 = 1.1*T;
A_new(id) = 0;
D_new = D;
for i=1:N
    % Double verification on artifacts or spikes
%Weighted Thresholding
if (i == 3 || i == 4 || i == 5 || i == 6) % D4, D5, D6 contains spike data, so high threshold
    k2 = 3;
else
    k2 = 1; % others more likely to be artifacts, so low threshold value
end

%%----verification ends here-------
sigma_sq = median(abs(D(i,:)))/0.6745;
% Th = (k2 - 0.1*i)*sqrt(2*log10(length(D))*sigma_sq^2);
Th(i) = k2*sqrt(2*log10(length(D))*sigma_sq^2);
idx = find((abs(D(i,:))> Th(i))==1);
tau = 0.01;
lamda1 = 1.1*Th(i);
%% Please uncomment below to include double verification of identifying artifact
if(abs(x_bpf(idx)) < 2*thr)  %| abs(x_hpf(idx)) < 2*thrh)
if(abs(x_bpf_sp(idx)) > thrp)
idx = [];
else
idx = idx;
end
else
idx = idx;
end
% if (i == 4 || i == 5 || i == 6)
% D_new(i,idx) = D(i,idx);
% else
% D_new(i,idx) = 0; % Hard
% D_new(i,idx) = sign(D(i,idx)).* abs(D(i,idx) - Th); % Soft
D_new(i,idx) = Th(i).^2./D(i,idx); % Garrote
% D_new(i,idx) = (sign(D(i,idx)).*(abs(D(i,idx)-Th)))./(1 + exp(- tau*(abs(D(i,idx)-lamda1)))); % SBSS
end
% end
%% Reconstruction
%%-SWT based reconstruction, ISWT
X_new = iswt(A_new, D_new, Lo_R, Hi_R); %X = ISWT(SWA(end,:),SWD,Lo_R,Hi_R)
data_new = X_new - mean(X_new);
%% Plot results
% plot(data_new);hold on;plot(sig_us,'r');
% xlabel('Time sample'); ylabel('Amplitude');
% legend('Data New', 'Artifactual Data');
id_A = id;
id_D = idx;
m = max(abs(A_old))/std(A_old);
Time_Required_Per_Sec_Data = toc/(data_length/Fs)

plot(sig_us,'k'); hold on;
plot(data_new,'r'); 
legend('Artifactual data', 'Cleaned data');
