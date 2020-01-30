%% Clearing the workspace and starting the timer

tic; % start the timer

clear all; close all; clc;


%% Loading the data

file = 'E:\LFPs\AEP2019\Rat62\22kSingleCall\';

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
             
numTrials = length(ts_events) % number of events

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
figure;
plot(tx, smooth(mean(lfp, 2), 35, 'loess'), '-k', 'LineWidth', 2)
xlim([-20 100])

%% Baseline period specification

% Define baseline period
baselinetime = [ -1000 00 ]; % in ms

% Convert baseline window time to indices
[~, baselineidx(1)] = min(abs(tx - baselinetime(1)));
[~, baselineidx(2)] = min(abs(tx - baselinetime(2)));

%% Parameters for wavelet transform

% Wavelet parameters
min_freq = 1;
max_freq = 150;
num_frex = 150;


% Other wavelet parameters
frex = linspace(min_freq, max_freq, num_frex);
time = -0.5:1/Fs:0.5;
half_wave = (length(time) - 1)/2;

% FFT parameters
nKern = length(time);
nData = size(lfp, 1)* size(lfp, 2);
nConv = nKern + nData - 1;

% Initialize output time-frequency data
tf = zeros(length(frex),size(lfp, 1));

%% Time-frequency decomposition

% FFT of total data
fft_lfp = fft( reshape(lfp, 1, []), nConv);


range_cycles = [ 2  10 ];
% cycles = logspace(log10(range_cycles(1)), log10(range_cycles(end)), num_frex);
cycles = linspace(range_cycles(1), range_cycles(end), num_frex);
num_cycles = length(range_cycles);

for cyclei = 1:length(num_cycles)
    
    for fi = 1:length(frex)
    
    % Create wavelet and get its FFT
    s        = cycles(fi) /(2*pi*frex(fi));

    wavelet  = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));
    
    % Run convolution for each of total, induced, and evoked
   

        % Need separate FFT 
        waveletX = fft(wavelet, nConv);
        waveletX = waveletX./max(waveletX);
        
        
        %  Notice that the fft_lfp cell changes on each iteration
        as = ifft(waveletX.*fft_lfp, nConv);
        as = as(half_wave + 1:end - half_wave);
        as = reshape(as, size(lfp, 1), size(lfp, 2));
        
            % Compute power
             tf(fi, :) = mean(abs(as).^2, 2);
             tf2(fi, :) = median(abs(as).^2, 2);
      
     % end loop around total, evoked, induced
end % end frequency loop

end % end cycle loop

   
%% Baseline corrected power

% dB converted power
baseline_power = mean(tf2(:, baselineidx(1):baselineidx(2)), 2);
dbconverted = 10*log10( bsxfun(@rdivide, tf2, baseline_power));

% Percent change in power w.r.t. baseline
pctchange = 100 * (tf2 - repmat(baseline_power, 1, size(lfp, 1)))./ repmat(baseline_power, 1, size(lfp, 1));

% Baseline division
baselinediv = tf2 ./ repmat(baseline_power, 1, size(lfp, 1));

% Z-transform
baseline_powerZ = tf2(:, baselineidx(1):baselineidx(2));
baselineZ = (tf2 - repmat(mean(baseline_powerZ, 2), 1, size(tf2, 2))) ./ repmat(std(baseline_powerZ, [], 2), 1, size(tf2, 2));

%% Plotting the results

figure('Color', [1 1 1]), clf
subplot(221)
contourf(tx, frex, dbconverted, 80, 'LineColor', 'none');
pcolor(tx, frex, dbconverted); shading interp;
colorbar;
xlabel('Time (ms)', 'FontSize', 14); ylabel('Frequency (Hz)', 'FontSize', 14);
title('dB converted and baseline corrected power')
set(gca, 'YDir', 'normal', 'XLim', [(x_min*Fs) + 100 (x_max*Fs) - 100], 'YLim', [min_freq max_freq], 'CLim', [0 2])
hold on
plot([0 0], [min_freq max_freq], '--k', 'LineWidth', 1.5)
hold off

subplot(222)
% contourf(tx, frex, pctchange, 80, 'LineColor', 'none');
pcolor(tx, frex, pctchange); shading interp;
colorbar;
xlabel('Time (ms)', 'FontSize', 14); ylabel('Frequency (Hz)', 'FontSize', 14);
title('% change in power from baseline')
set(gca, 'YDir', 'normal', 'XLim', [(x_min*Fs) + 100 (x_max*Fs) - 100], 'YLim', [min_freq max_freq], 'CLim', [-50 50])
hold on
plot([0 0], [min_freq max_freq], '--k', 'LineWidth', 1.5)
hold off

subplot(223)
% contourf(tx, frex, dbconverted, 80, 'LineColor', 'none'); 
pcolor(tx, frex, dbconverted); shading interp;
colorbar;
xlabel('Time (ms)', 'FontSize', 14); ylabel('Frequency (Hz)', 'FontSize', 14);
title('baseline divided power')
set(gca, 'YDir', 'normal', 'XLim', [(x_min*Fs) + 100 (x_max*Fs) - 100], 'YLim', [min_freq max_freq], 'CLim', [0 2])
hold on
plot([0 0], [min_freq max_freq], '--k', 'LineWidth', 1.5)
hold off

subplot(224)
% contourf(tx, frex, dbconverted, 80, 'LineColor', 'none'); 
pcolor(tx, frex, dbconverted); shading interp;
colorbar;
xlabel('Time (ms)', 'FontSize', 14); ylabel('Frequency (Hz)', 'FontSize', 14);
title('z transformed power')
set(gca, 'YDir', 'normal', 'XLim', [(x_min*Fs) + 100 (x_max*Fs) - 100], 'YLim', [min_freq max_freq], 'CLim', [0 2])
hold on
plot([0 0], [min_freq max_freq], '--k', 'LineWidth', 1.5)
hold off

%% Mean vs median

% db-correct and plot
labelz = {'mean';'median'};


    mean_baseline_power = mean(tf(:, baselineidx(1):baselineidx(2)), 2);
    mean_dbconverted = 10*log10( bsxfun(@rdivide, tf(:, :), mean_baseline_power));
    median_baseline_power = mean(tf2(:, baselineidx(1):baselineidx(2)), 2);
    median_dbconverted = 10*log10( bsxfun(@rdivide, tf2(:, :), median_baseline_power));
    
    % Plot
    figure
    subplot(2, 2, 1)
%     contourf(tx, frex, mean_dbconverted, 80, 'LineColor', 'none')
    pcolor(tx, frex, mean_dbconverted); shading interp;
    set(gca, 'CLim',[0 2], 'XLim', [(x_min*Fs) + 100 (x_max*Fs) - 100], 'YLim', [min_freq max_freq], 'YScale', 'linear')
    title(labelz{1}), ylabel('Frequency (Hz)'), xlabel('Time (ms)')
    hold on
    plot([0 0], [0 90], '--k', 'LineWidth', 1.5)
    hold off

    subplot(2,2,2)
%     contourf(tx,frex, median_dbconverted, 80, 'Linecolor', 'none')
    pcolor(tx, frex, median_dbconverted); shading interp;
    set(gca, 'CLim',[0 2], 'XLim',[(x_min*Fs) + 100 (x_max*Fs) - 100], 'YLim', [min_freq max_freq], 'YScale', 'linear')
    title(labelz{2}), ylabel('Frequency (Hz)'), xlabel('Time (ms)')
    hold on
    plot([0 0], [0 max_freq], '--k', 'LineWidth', 1.5)
    hold off

% Plot relationship between mean and median
subplot(223)
db_mean = 10*log10( bsxfun(@rdivide,tf(:, :), mean(tf(:, baselineidx(1):baselineidx(2)), 2)));
db_medn = 10*log10( bsxfun(@rdivide,tf2(:, :), mean(tf2(:, baselineidx(1):baselineidx(2)), 2)));
plot(db_mean(:),db_medn(:), '.')
r = corr(db_mean(:), db_medn(:));
legend([ 'R^2 = ' num2str(r*r) ])
xlabel('dB from Mean')
ylabel('dB from Median')

%% Signal-to-noise ratio calculation

snr_bs = zeros(length(frex), size(lfp, 1));
snr_tf = zeros(length(frex), size(lfp, 1));
tf     = zeros(length(frex), size(lfp, 1));
wavelet_cycles = 4;

for fi = 1:length(frex)
    
    % create wavelet and get its FFT
    wavelet = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*frex(fi)))^2))/frex(fi);
    fft_wavelet = fft(wavelet, nConv);
    
    % run convolution
    convolution_result = ifft(fft_wavelet.*fft_lfp, nConv) * sqrt(wavelet_cycles /(2*pi*frex(fi)));
%     convolution_result = convolution_result(1:nConv);
    convolution_result = convolution_result(half_wave + 1:end - half_wave);
    convolution_result = reshape(convolution_result, size(lfp, 1), size(lfp, 2));
    
    % extract SNR in two ways
    snr_tf(fi, :) = mean(abs(convolution_result).^2, 2)./std(abs(convolution_result).^2, [], 2);
    snr_bs(fi, :) = mean(abs(convolution_result).^2, 2)./std(mean(abs(convolution_result(baselineidx(1):baselineidx(2), :)).^2, 1), [], 2);
    
    % and extract trial-averaged power
    tf(fi, :) = mean(abs(convolution_result).^2, 2);
    
end

% Plot
figure
% subplot(121)
% contourf(tx, frex, snr_bs, 180, 'LineColor', 'none')
% set(gca, 'CLim', [0.5 2], 'XLim', [-200 1000], 'YLim', [frex(1) max_freq])
% colorbar
% title('SNR_b_a_s_e_l_i_n_e (mean/std)'), ylabel('Frequency (Hz)'), xlabel('Time (ms)')
% axis square


subplot(122)


baseline_power = mean(tf(:,baselineidx(1):baselineidx(2)),2);
baselinediv    = tf ./ repmat(baseline_power, 1, size(lfp, 1));

plot(snr_bs(1:3:end), baselinediv(1:3:end), '.')
xlabel('SNR_b_a_s_e_l_i_n_e'), ylabel('Power (/baseline)')
axis square


% figure
% subplot(121)
% contourf(tx, frex, snr_tf, 80, 'LineColor', 'none')
% set(gca,'clim',[.5 1.25],'xlim',[-200 1000],'ylim',[frex(1) max_freq])

% colorbar
title('SNR_t_f (mean/std)'), ylabel('Frequency (Hz)'), xlabel('Time (ms)')
axis square

subplot(122)
plot(snr_tf(1:3:end), baselinediv(1:3:end), '.')
xlabel('SNR_t_f'), ylabel('Power (/baseline)')
axis square

% Time-series of SNR
figure
plot(tx, smooth(abs(mean(lfp(:, :), 2) ./ std(lfp(:, :), [], 2)), 50, 'loess'))
title('Time-domain SNR time series')
set(gca,'XLim',[(x_min*Fs) + 100 (x_max*Fs) - 100])
xlabel('Time (ms)'), ylabel('SNR')

% Now compute SNR of peak compared to prestim noise
stimeidx = dsearchn(tx', 150);
etimeidx = dsearchn(tx', 400);
% disp([ 'ERP SNR between 150 and 400 ms at FCz: ' num2str(max(mean(lfp(stimeidx:etimeidx,:), 2)) / std(mean(lfp(baselineidx(1):baselineidx(2), :), 2),[], 2)) ])

% save(strcat(file, 'lfpBLA_04_40cycles_1000ms'), 'lfp', 'tx')
% save(strcat(file, 'SpectralPowerBLA_04_40cycles_1000ms'), 'baselinetime', 'dbconverted', 'tx', 'frex', 'tf', 'tf2')
toc;

%% end of code