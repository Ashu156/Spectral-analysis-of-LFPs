clear all; close all; clc;

%% Loading the data

file = 'E:\LFPs\AEP2019\Rat59\22kSingleCall\';

load(strcat(file, 'matfile.mat'));

x_min = -1.5; x_max = 2.5;

% Choosing the required channel

ts_csc = CSC28_TS; % Channel of choice
dp_csc = CSC28_DP; % Data points in channel of choic0
ADBitVolts = (str2double(CSC28_NlxHeader{15}(13:38)))*10^6; % in microVolts
Fs = round(str2double(CSC28_NlxHeader{14}(end - 4:end))); % sampling frequency

if Fs > 1e3
    Fs = 1e3;
end

dp_csc = dp_csc*ADBitVolts; % CSC data points in microVolts
dp_csc = dp_csc(:);

dp_csc = detrend(dp_csc,'constant');
ts_events = (Events_tone_on - 0000.0)./10^6; % selecting events of choice (onset of white noise pips in this case)

% ts_events = [ts_events(18:31) ts_events(33:48) ts_events(55:57) ts_events(59) ts_events(61:62)...
%              ts_events(64:67) ts_events(69:70) ts_events(72:73) ts_events(76:78) ts_events(80:88) ts_events(90)...                         
%              ts_events(93:94) ts_events(97:98) ts_events(114) ts_events(116:129) ts_events(131:134)...
%              ts_events(136:137) ts_events(140) ts_events(143:145) ts_events(147:151) ts_events(156:166)...
%              ts_events(168)... % ts_events(143:145) ts_events(147:151) ts_events(156:166)...
%              ts_events(168) ts_events(170) ts_events(172) ts_events(174:178) ts_events(180:190)...
%              ts_events(204:212) ts_events(214) ts_events(216) ts_events(218) ts_events(220)...
%              ts_events(222:227) ts_events(229:233) ts_events(235) ts_events(237:241) ts_events(244:254)...
%              ts_events(256:268) ts_events(270:275)];
             
numTrials = length(ts_events); % number of events

% if numTrials > 100
%     numTrials = 100;
% end

ts_csc = ts_csc./10^6 ; % in seconds

% Generate all the timestamps corresponding to all data points
tt = [0:1/Fs:511/Fs]' ;
tts_csc = [];
for i = 1:length(ts_csc)
    tts_csc = [tts_csc  [ts_csc(i) + tt]];
end
tts_csc = tts_csc(:);
kmin = round(Fs*(x_min));         % Min x-limit  
kmax = round(Fs*(x_max));         % Max x-limit  

tx = (x_min:1/Fs:x_max)*1000 ; % time in milliseconds


lfp = zeros(length([kmin:kmax]),numTrials);

for i = 1:numTrials 
    min_dist(1, i) = min(abs(tts_csc(:) - ts_events(i)));
    k0(1, i) = find(abs(tts_csc(:) - ts_events(i)) == min_dist(1,i)) ;
    lfp(:, i) = dp_csc(k0(1, i) + kmin : k0(1, i) + kmax);
end

lfp = lfp(1:length(tx), :);

%%

baseline = [-100 0];
baselineidx = dsearchn(tx', baseline');

for i = 1:numTrials
    
    baselinedata(:, i) = lfp(baselineidx(1):baselineidx(2), i);
    
end

ashu = repmat(baselinedata, 4, 1);
ashuMean = mean(ashu, 2);

for i = 1:numTrials
    erp(:, i) = lfp(baselineidx(1):baselineidx(1) + size(ashu, 1) , i);
end

erp = erp(1:end - 1, :);
aep = erp - ashu;

for i = 1:numTrials
    aep2(:, i) = erp(:, i) - ashuMean;
end

figure;
plot(smooth(median(aep, 2), 15*2, 'sgolay'), 'linew', 2)
% hold on
% plot([100 100], [-400 600], '--k', 'linew', 1.5);
xlim([0 200])

figure; 
plot(smooth(median(aep2, 2), 15*2, 'sgolay'), 'linew', 2)
% hold on
% plot([100 100], [-400 600], '--k', 'linew', 1.5);
xlim([0 200])

%% end of code
