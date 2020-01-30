%% Clearing the workspace and starting the timer

tic; % start the timer

clear all; close all; clc;


%% Loading the data

file = 'E:\LFPs\AEP2019\Rat81\10kpips\';

load(strcat(file, 'matfile.mat'));

x_min = -0.5; x_max = 2.5;

% Choosing the required channel

ts_csc = CSC02_TS; % Channel of choice
dp_csc = CSC02_DP; % Data points in channel of choic0
ADBitVolts = (str2double(CSC02_NlxHeader{15}(13:38)))*10^6; % in microVolts
Fs = round(str2double(CSC02_NlxHeader{14}(end - 4:end))); % sampling frequency

if Fs > 1e3
    Fs = 1e3;
end

dp_csc = dp_csc*ADBitVolts; % CSC data points in microVolts
dp_csc = dp_csc(:);

dp_csc = detrend(dp_csc,'constant');
ts_events = (Events_tone_on - 0000.0)./10^6; % selecting events of choice (onset of white noise pips in this case)

% ts_events = [ts_events(8:37) ts_events(43:48) ts_events(50:56) ts_events(58:68) ts_events(70:73)...
%              ts_events(75:76) ts_events(89:91) ts_events(94:98) ts_events(100:108) ts_events(110:122) ts_events(125:130)...                         
%              ts_events(132:133) ts_events(135:150) ts_events(152) ts_events(154:158) ts_events(175:177)...
%              ts_events(179:180) ts_events(182:185) ts_events(187:190) ts_events(192:208) ts_events(215:239)...
%              ts_events(241:245)];% ts_events(143:145) ts_events(147:151) ts_events(156:166)...
%              ts_events(168) ts_events(170) ts_events(172) ts_events(174:178) ts_events(180:190)...
%              ts_events(204:212) ts_events(214) ts_events(216) ts_events(218) ts_events(220)...
%              ts_events(222:227) ts_events(229:233) ts_events(235) ts_events(237:241) ts_events(244:254)...
%              ts_events(256:268) ts_events(270:275)];
             
numTrials = length(ts_events) % number of events

if numTrials > 100
    numTrials = 100
end

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
save(strcat(file, 'lfpBLA_26_10cycles_1000ms'), 'lfp', 'tx')
figure;
plot(tx, smooth(mean(lfp, 2), 35, 'loess'), '-k', 'LineWidth', 2)
xlim([-20 100])

%% Baseline period specification

% Define baseline period
baselinetime = [ -100 00 ]; % in ms

% Convert baseline window time to indices
[~, baselineidx(1)] = min(abs(tx - baselinetime(1)));
[~, baselineidx(2)] = min(abs(tx - baselinetime(2)));

%%

low_freq = 0;
high_freq = 150;
fs = 1000;  % Sample rate
xtime  = tx;
tmin = xtime(1); tmax = xtime(end);
np = length(xtime);
freq = max([low_freq 2]):high_freq;  % Every Hz, starting from at least 1Hz

%% Time-frequency decomposition

for i = 1:numTrials
    
%    TFmapmap = qcwt( lfp(:, i), fs, freq,  5/(2*sqrt(2*log(2))));
%    TF = qcwt( lfp(:, i), fs, freq,  5/(2*sqrt(2*log(2))));
   TF = qcwt( lfp(:, i), fs, freq,  2);
%    TF = qcwt( lfp(:, i), fs, freq,  [2:0.1:10], 'Method', 'MinIP');
   TFmap(i, :, :) = TF;
   
end


   
%% Baseline corrected power

ashu = squeeze(mean(TFmap, 1));
% dB converted power
baseline_power = mean(ashu(:, baselineidx(1):baselineidx(2)), 2);
dbconverted = 10*log10( bsxfun(@rdivide, ashu, baseline_power));

% Percent change in power w.r.t. baseline
pctchange = 100 * (ashu - repmat(baseline_power, 1, size(lfp, 1)))./ repmat(baseline_power, 1, size(lfp, 1));

% Baseline division
baselinediv = ashu ./ repmat(baseline_power, 1, size(lfp, 1));

% Z-transform
baseline_powerZ = ashu(:, baselineidx(1):baselineidx(2));
baselineZ = (ashu - repmat(mean(baseline_powerZ, 2), 1, size(ashu, 2))) ./ repmat(std(baseline_powerZ, [], 2), 1, size(ashu, 2));

%% Plotting the results

figure('Color', [1 1 1]), clf
pcolor(xtime, freq, dbconverted);
shading interp
colormap jet;
colorbar;
xlabel('Time (ms)', 'FontSize', 14); ylabel('Frequency (Hz)', 'FontSize', 14);
title('dB converted and baseline corrected power')
% set(gca, 'YDir', 'normal', 'XLim', [(x_min*Fs) + 100 (x_max*Fs) - 100], 'YLim', [min_freq max_freq], 'CLim', [0 2])
hold on
plot([0 0], [1 80], '--w', 'LineWidth', 1.5)
hold off
xlim([-200 1000]); ylim([1 80])

figure('Color', [1 1 1]), clf
pcolor(xtime, freq, pctchange); colorbar;
shading interp; colormap jet;
xlabel('Time (ms)', 'FontSize', 14); ylabel('Frequency (Hz)', 'FontSize', 14);
title('% change in power from baseline')

hold on
plot([0 0], [1 80], '--w', 'LineWidth', 1.5)
hold off
xlim([-200 1000]); ylim([1 80])

figure('Color', [1 1 1]), clf
pcolor(xtime, freq, baselinediv); colorbar;
shading interp; colormap jet; 
xlabel('Time (ms)', 'FontSize', 14); ylabel('Frequency (Hz)', 'FontSize', 14);
title('baseline divided power')
% set(gca, 'YDir', 'normal', 'XLim', [(x_min*Fs) + 100 (x_max*Fs) - 100], 'YLim', [min_freq max_freq], 'CLim', [0 2])
hold on
plot([0 0], [1 80], '--w', 'LineWidth', 1.5)
hold off
xlim([-200 1000]); ylim([1 80])

figure('Color', [1 1 1]), clf
pcolor(xtime, freq, baselineZ); colorbar;
shading interp; colormap jet;
xlabel('Time (ms)', 'FontSize', 14); ylabel('Frequency (Hz)', 'FontSize', 14);
title('z transformed power')
% set(gca, 'YDir', 'normal', 'XLim', [(x_min*Fs) + 100 (x_max*Fs) - 100], 'YLim', [min_freq max_freq], 'CLim', [0 2])
hold on
plot([0 0], [1 80], '--k', 'LineWidth', 1.5)
hold off
xlim([-200 1000]); ylim([1 80])


toc;

%% Plotting surf and its projection simultaneously

min_x = min(min(xtime));
min_y = min(min(freq));
max_x = max(max(xtime));
max_y = max(max(freq));

figure; hold on;
surf(xtime, freq, baselineZ);
shading interp

planeimg = abs(baselineZ);

imgzposition = -40;

surf([min_x max_x],[min_y max_y],repmat(imgzposition, [2 2]),...
planeimg,'facecolor','texture')
colormap(jet);
view([-38 30])
caxis([-10 20]); colorbar;
%% end of code