%% Description:
% 
% This script loads all the analyzed data in the workspace and further plots
% the AEPs and power spectral density for induvidual rats averaged across trials. 
% Further, AEPs and power spectral density averaged across all the rats and
% trials is plotted in a separate figure.
% 
%  The script is divided in two halves: the first half deals with LFPs
%  recorded from BLA, ehile the later half deal swith LFPs recorded from
%  dmPFC. Each half is further divided into two sub-halves. While the first
%  sub-half deals with LFPs recorded while a single 22-kHz USV call was
%  presented to the rat, the second half deals with LFPs recorded while
%  rats were subjected to 10-kHz pips. 
% 
% N.B.: One should be careful about which cells in the script are uncommented!!!

%% 0.

clear all; close all; clc;

%% 1. Rat IDs for 22kSingleCall BLA

fileLocation = 'E:\LFPs\AEP2019\Rat';

% Rat IDs for 22kSingleCall

 ratID        = {'54'; '57'; '58';'59'; '62'; ...
                 '64'; '67'; '75';'81'; '82'}; %             
                  
 exptID       = '22kSingleCall';
 fileID       = 'lfp';
 brainArea    = 'BLA';
 nCycles      = '40cycles'; 
 
% Selected channels for 22kSingleCall
               
  channelID    = { '28'; '28'; '26'; '28'; '28'; ...
                  '28'; '28'; '28'; '01'; '29'};
 
 

for i = 1:length(ratID)   % 
    load(strcat(fileLocation, ratID{i, 1}, '\', exptID, '\', fileID, brainArea, '_',  channelID{i, 1}, '_', nCycles, '_1000ms', '.mat'));
% load(strcat(fileLocation, ratID{i, 1}, '\', exptID, '\', fileID, brainArea, '_',  channelID{i, 1},  '_', nCycles, '_100ms', '.mat'));
end

% clear fileLocation ratID exptID fileID brainArea channelID 

%% 2. Rat IDs for 10kpips BLA

fileLocation = 'E:\LFPs\AEP2019\Rat';

% Rat IDs for 10kpips            
 ratID        = {'54'; '58'; '59'; '62';... 
                 '64';'75'; '81'; '82'; };       
             
                  
 exptID       = '10kpips';
 fileID       = 'lfp';
 brainArea    = 'BLA';
 nCycles      = '40cycles';

 
%  Selected channels for 10kpips
channelID    = {'28'; '26'; '28'; '28';...
                '28';'28'; '01'; '29'; };
            
channelID    = {'29'; '28'; '29'; '29';...
                '29';'29'; '02'; '31';}; 

for i = 1:length(ratID)
    load(strcat(fileLocation, ratID{i, 1}, '\', exptID, '\', fileID, brainArea, '_',  channelID{i, 1}, '_', nCycles, '_100ms', '.mat'));
end

% clear fileLocation ratID exptID fileID brainArea channelID i

%% 3. Mean AEP in response to 22 kHz single call / 10 kHz pips

% %%%%%%%%%%%%%%%%%%%%%%%% -----  BLA ----- %%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(exptID,'22kSingleCall') == 1
    
    lfp = zeros(4001, 100);
    lfp_z = zeros(4001, 100);
    
else
    
    lfp = zeros(3001, 100);
    lfp_z = zeros(3001, 100);
    
end


% Determining the number of rows and columns in the figure
nCol = 5;
Q = floor(length(ratID)/nCol);
R = mod(length(ratID), nCol);

if R == 0
    nRow = Q;
else
    nRow = Q + 1;
end
    
for i = 1:length(ratID) % 1
    
    lfp = evalin('base', strcat('lfp_', ratID{i}));
    lfp_z = (zscore(lfp, 0, 1));
    subplot(nRow, nCol, i)
    plot(tx_54, smooth(mean(lfp_z, 2), 35, 'sgolay'), '-k', 'LineWidth', 2);
    hold on
    xlim([-20 100])
%     plot([0 0], [-0.1 0.2], '--r', 'LineWidth', 2)
    set(gcf, 'Color', [1 1 1])
    ylabel('AEP amplitude (z-score)')
    xlabel('Time (ms)')
    title(strcat('Rat', ratID{i}))
    
end

figure;

lfpAll = zeros(length(ratID), 4001, 100);
lfpzAll = zeros(length(ratID), 4001, 100);

for i = 1:length(ratID) 
    lfp = evalin('base', strcat('lfp_', ratID{i}));
    lfpAll(i, :, :) = lfp;
    lfp_z = (zscore(lfp(400:600, :) , 0, 1));
    lfpzAll = lfp_z;
    plot(-100:100, smooth(mean(lfp_z, 2), 35, 'sgolay'), 'LineWidth', 2);
    hold on
    xlim([-20 100])
    plot([0 0], [-0.2 0.2], '--k', 'LineWidth', 1.5)
    set(gcf, 'Color', [1 1 1])
    ylabel('AEP amplitude (z-score)')
    xlabel('Time (ms)')
    title('AEPs from each rat averaged across all trials')
end

lfp = lfp./length(ratID);
lfp_z = lfp_z./length(ratID);


   
%% 4. Mean power spectral density in response to 22 kHz single call

% %%%%%%%%%%%%%%%%%%%%%%%% -----  BLA ----- %%%%%%%%%%%%%%%%%%%%%%%%


 fileLocation = 'E:\LFPs\AEP2019\Rat';
 
 % Rat IDs for 22kSingleCall 
 
 ratID        = {'54'; '57'; '58';'59'; '62'; ...
                 '64'; '67'; '75';'81'; '82'}; %  
                 
            
                  
 exptID       = '22kSingleCall';
 fileID       = 'SpectralPower';
 brainArea    = 'BLA';
 nCycles      = '40cycles';
 
 % Selected channels for 22kSingleCall

 channelID    = { '28'; '28'; '26'; '28'; '28'; ...
                  '28'; '28'; '28'; '01'; '29'};

 
 
for i = 1:length(ratID)
    load(strcat(fileLocation, ratID{i, 1}, '\', exptID, '\', fileID, brainArea, '_',  channelID{i, 1}, '_',  nCycles, '_1000ms', '.mat'));
end

tf2all = zeros(150, 4001);


% Determining the number of rows and columns in the figure

% figure;
nCol = 5;
Q = floor(length(ratID)/nCol); 
R = mod(length(ratID), nCol);

if R == 0
    nRow = Q;
else
    nRow = Q + 1;
end

% Baseline for normalization
baselinetime = [ -100 0 ];
baselineidx = dsearchn(tx_54', baselinetime');
dbconverted = [];


for i = 1:length(ratID)    % 
    
    tf2 = evalin('base', strcat('tf2_', ratID{i}));
    
    % dB converted power
    baseline_power = mean(tf2(:, baselineidx(1):baselineidx(2)), 2);
    dbconverted{i, 1} = 10*log10( bsxfun(@rdivide, tf2, baseline_power));
    
    % Percent change in power w.r.t. baseline
    pctchange{i, 1} = 100 * (tf2 - repmat(baseline_power, 1, size(lfp, 1)))./ repmat(baseline_power, 1, size(lfp, 1));
    
    % Baseline division
    baselinediv{i, 1} = tf2 ./ repmat(baseline_power, 1, size(lfp, 1));

    % Z-transform
    baseline_powerZ = tf2(:, baselineidx(1):baselineidx(2));
    baselineZ{i, 1} = (tf2 - repmat(mean(baseline_powerZ, 2), 1, size(tf2, 2))) ./ repmat(std(baseline_powerZ, [], 2), 1, size(tf2, 2));

%     tf2all = tf2all + tf2;
%     tf2all = tf2all + dbconverted{i, 1};
%     tf2all = tf2all + pctchange{i, 1};
%     tf2all = tf2all + baselinediv{i, 1};
    tf2all = tf2all + baselineZ{i, 1};
%     subplot(nRow, nCol, i)
%     imagesc(tx_54, frex_54, dbconverted); colorbar;
%    figure(i)
%     contourf(tx_54, frex_54, baselineZ{i, 1}, 200, 'LineColor', 'none'); colorbar;
%     hold on
%     plot([0 0], [1 150], '--k', 'LineWidth', 3);
%     plot([1300 1300], [1 150], '--k', 'LineWidth', 3);
%     set(gca, 'YDir', 'normal', 'CLim', [0 5]);
%     set(gcf, 'Color', [1 1 1]);
%     colormap parula
%     xlabel('Time (ms)')
%     ylabel('Frequency (Hz)')
%     set(gca, 'XLim', [-1000 2400],  'XTick', [-1000:500:2400], 'CLim', [0 5])
    
%     title(strcat('Rat', ratID{i}, '-', nCycles, '-', 'CSC', channelID{i}))
    
end

tf2all = tf2all./length(ratID);
figure;
% contourf(tx_54, frex_54, tf2all, 200, 'LineColor', 'none');
pcolor(tx_54, frex_54, tf2all); shading interp;
hold on
plot([0 0], [1 150], '--k', 'LineWidth', 1.5);
set(gca, 'CLim', [-1.5 1.5])
% colormap jet
xlabel('Time(ms)'); ylabel('Frequency (Hz)');
xlim([-1000 2400])
title('Power spectral density averaged across all rats and trials')

%% Contour plots for individual rats

tf2 = evalin('base', strcat('tf2_', ratID{5}));
    % dB converted power
    baseline_power = mean(tf2(:, baselineidx(1):baselineidx(2)), 2);
    dbconverted_smpl = 10*log10( bsxfun(@rdivide, tf2, baseline_power));
    
    % Percent change in power w.r.t. baseline
    pctchange_smpl = 100 * (tf2 - repmat(baseline_power, 1, size(lfp, 1)))./ repmat(baseline_power, 1, size(lfp, 1));
    
    % Baseline division
    baselinediv_smpl = tf2 ./ repmat(baseline_power, 1, size(lfp, 1));

    % Z-transform
    baseline_powerZ = tf2(:, baselineidx(1):baselineidx(2));
    baselineZ_smpl = (tf2 - repmat(mean(baseline_powerZ, 2), 1, size(tf2, 2))) ./ repmat(std(baseline_powerZ, [], 2), 1, size(tf2, 2));

    figure(1)
%     contourf(tx_54, frex_54, baselineZ_smpl, 200, 'LineColor', 'none')
    pcolor(tx_54, frex_54, baselineZ_smpl); shading interp
    set(gca, 'YDir', 'normal')
    caxis([-5 5])
    xlim([-1000 2400])
    ylim([1 50])
    
%% 5. Mean power spectral density in response to 10 kHz pips

% %%%%%%%%%%%%%%%%%%%%%%%% -----  BLA ----- %%%%%%%%%%%%%%%%%%%%%%%%


 fileLocation = 'E:\LFPs\AEP2019\Rat';
 
 
 % Rat IDs for 10kpips    
 
 ratID        = {'54'; '58'; '59'; '62'; '64'; '71';...
                 '75'; '81'; '82'; '84'; '86'; '88'}; 
                        
                             
 exptID       = '10kpips';
 fileID       = 'SpectralPower';
 brainArea    = 'BLA';
 nCycles      = '40cycles';
 
 %  Selected channels for 10kpips
channelID    = {'28'; '26'; '28'; '28'; '28'; '28';...
                '28'; '01'; '29'; '29'; '26'; '28'};
            
channelID    = {'29'; '28'; '29'; '29'; '29'; '29';...
                '29'; '02'; '31'; '31'; '28'; '29'}; 

for i = 1:length(ratID)
    load(strcat(fileLocation, ratID{i, 1}, '\', exptID, '\', fileID, brainArea, '_',  channelID{i, 1}, '_', nCycles, '_100ms', '.mat'));
end
tf2all = zeros(150, 3001);


% Determining the number of rows and columns in the figure
figure;
nCol = 5;
Q = floor(length(ratID)/nCol); 
R = mod(length(ratID), nCol);

if R == 0
    nRow = Q;
else
    nRow = Q + 1;
end

% Baseline for normalization
baselinetime = [-100 0];
baselineidx = dsearchn(tx_54', baselinetime');

for i = 1:length(ratID)   % 
    
    tf2 = evalin('base', strcat('tf2_', ratID{i}));
    
    % dB converted power
    baseline_power = mean(tf2(:, baselineidx(1):baselineidx(2)), 2);
    dbconverted{i, 1} = 10*log10( bsxfun(@rdivide, tf2, baseline_power));
    
    % Percent change in power w.r.t. baseline
    pctchange{i, 1} = 100 * (tf2 - repmat(baseline_power, 1, size(lfp, 1)))./ repmat(baseline_power, 1, size(lfp, 1));
    
    % Baseline division
    baselinediv{i, 1} = tf2 ./ repmat(baseline_power, 1, size(lfp, 1));

    % Z-transform
    baseline_powerZ = tf2(:, baselineidx(1):baselineidx(2));
    baselineZ{i, 1} = (tf2 - repmat(mean(baseline_powerZ, 2), 1, size(tf2, 2))) ./ repmat(std(baseline_powerZ, [], 2), 1, size(tf2, 2));

%     tf2all = tf2all + tf2;
%     tf2all = tf2all + dbconverted{i, 1};
%     tf2all = tf2all + pctchange{i, 1};
%     tf2all = tf2all + baselinediv{i, 1};
    tf2all = tf2all + baselineZ{i, 1};
%     subplot(nRow, nCol, i)
%     imagesc(tx_54, frex_54, dbconverted); colorbar;
%    figure(i)
%     contourf(tx_54, frex_54, baselineZ{i, 1}, 200, 'LineColor', 'none'); colorbar;
%     hold on
%     plot([0 0], [1 150], '--k', 'LineWidth', 3);
%     plot([1300 1300], [1 150], '--k', 'LineWidth', 3);
%     set(gca, 'YDir', 'normal', 'CLim', [0 5]);
%     set(gcf, 'Color', [1 1 1]);
%     colormap parula
%     xlabel('Time (ms)')
%     ylabel('Frequency (Hz)')
%     set(gca, 'XLim', [-1000 2400],  'XTick', [-1000:500:2400], 'CLim', [0 5])
    
%     title(strcat('Rat', ratID{i}, '-', nCycles, '-', 'CSC', channelID{i}))
    
end

tf2all = tf2all./length(ratID);
figure;
contourf(tx_54, frex_54, tf2all, 200, 'LineColor', 'none');
hold on
plot([0 0], [1 150], '--k', 'LineWidth', 1.5);
set(gca, 'CLim', [-1.5 1.5])
% colormap jet
xlabel('Time(ms)'); ylabel('Frequency (Hz)');
xlim([-1000 2400])
title('Power spectral density averaged across all rats and trials')

%% SECOND HALF OF SCRIPT STARTS HERE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ----- dmPFC----- %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. Rat IDs for 22kSingleCall dmPFC

fileLocation = 'E:\LFPs\AEP2019\Rat';


% Rat IDs for 22kSingleCall
 ratID        = {'26'; '28'; '29'; '30'; '32'; '37'; '39'; '45'; '47'; '48';... 
                 '49'; '50'; '51'; '52'; '53'; '54';'55'; '56'; '59'; '61';...                  
                 '62'; '63'; '64'; '66';'67'; '71'; '73'; '74'; '75'; '76'};
             
 exptID       = '22kSingleCall';
 fileID       = 'lfp';
 brainArea    = 'dmPFC';
 
 
% Selected channels for 22kSingleCall

 channelID    = {'31'; '31'; '31'; '31'; '31'; '30'; '30'; '30'; '30'; '29'; '30';...
                 '30'; '30'; '30'; '30'; '30'; '30'; '30'; '30'; '30'; '31';... 
                 '30'; '30'; '30'; '30'; '30'; '30'; '30'; '30'; '30'};
             
%  channelID    = {'32'; '32'; '32'; '32'; '32'; '31'; '31'; '31'; '31'; '30'; '31';...
%                  '31'; '31'; '31'; '31'; '31'; '31'; '31'; '31'; '31'; '32';... 
%                  '31'; '31'; '31'; '31'; '31'; '31'; '31'; '31'; '31'; '31'};
             
%  channelID    = {'31'; '31'; '31'; '31'; '31'; '30'; '30'; '30'; '30'; '29'; '30';...
%                  '30'; '30'; '30'; '30'; '30'; '30'; '30'; '30'; '30'; '31';... 
%                  '30'; '30'; '30'; '30'; '30'; '30'; '30'; '30'; '30'; '30'};
 
 

for i = 1:length(ratID)
    load(strcat(fileLocation, ratID{i, 1}, '\', exptID, '\', fileID, brainArea, '_',  channelID{i, 1}, '.mat'));
end

% clear fileLocation ratID exptID fileID brainArea channelID i

%% 2. Rat IDs for 10kpips dmPFC

fileLocation = 'E:\LFPs\AEP2019\Rat';


% Rat IDs for 10kpips

 ratID        = {'26'; '28'; '29'; '30'; '32'; '39'; '45'; '47'; '48';... 
                 '49'; '51'; '52'; '53'; '54';'55'; '56'; '59'; '61'; '62';...                  
                 '63'; '64';'71'; '73'; '74'; '75'; '76'; '77'};      
             
                  
 exptID       = '10kpips';
 fileID       = 'lfp';
 brainArea    = 'dmPFC';
 

 
%  Selected channels for 10kpips

%  channelID    = {'31'; '31'; '31'; '31'; '31'; '30'; '30'; '30'; '29'; '30';...
%                  '30'; '30'; '30'; '30'; '30'; '30'; '30'; '30'; '31';... 
%                  '30'; '30'; '30'; '30'; '30'; '30'; '30'};
             
%  channelID    = {'32'; '32'; '32'; '32'; '32'; '31'; '31'; '31'; '30'; '31';...
%                  '31'; '31'; '31'; '31'; '31'; '31'; '31'; '31'; '32';... 
%                  '31'; '31'; '31'; '31'; '31'; '31'; '31'};
% 
channelID    = {'31'; '31'; '31'; '31'; '31'; '30'; '30'; '30'; '29'; '31';...
                 '30'; '30'; '30'; '30'; '30'; '30'; '30'; '30'; '31';... 
                 '30'; '30'; '30'; '30'; '30'; '30'; '30'; '30'};

for i = 1:length(ratID)
    load(strcat(fileLocation, ratID{i, 1}, '\', exptID, '\', fileID, brainArea, '_',  channelID{i, 1}, '.mat'));
end

% clear fileLocation ratID exptID fileID brainArea channelID i


 %%  3. Mean AEP in response to 22 kHz single call / 10kHz pips

% %%%%%%%%%%%%%%%%%%%%%%%% ----- dmPFC ----- %%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(exptID,'22kSingleCall') == 1
    
    lfp = zeros(2001, 100);
    lfp_z = zeros(2001, 100);
    
else
    
    lfp = zeros(2001, 180);
    lfp_z = zeros(2001, 180);
    
end


% Determining the number of rows and columns in the figure
nCol = 5*2;
Q = floor(length(ratID)/nCol);
R = mod(length(ratID), nCol);

if R == 0
    nRow = Q;
else
    nRow = Q + 1;
end
    
for i = 1:length(ratID)
    
    lfp = lfp + evalin('base', strcat('lfp_', ratID{i}));
    lfp_z = (zscore(lfp, 0, 1));
    subplot(nRow, nCol, i)
    plot(tx_54, smooth(mean(lfp_z, 2), 30, 'sgolay'), '-k', 'LineWidth', 3);
    hold on
    xlim([-20 200])
    plot([0 0], [-0.1 0.1], '--r', 'LineWidth', 2)
    set(gcf, 'Color', [1 1 1])
    ylabel('z-scored AEP')
    xlabel('Time (ms)')
    title(strcat('Rat', ratID{i}))
    
end

figure;

for i = 1:length(ratID)
    lfp = lfp + evalin('base', strcat('lfp_', ratID{i}));
    lfp_z = (zscore(lfp, 0, 1));
    plot(tx_54, smooth(mean(lfp_z, 2), 30, 'sgolay'), 'LineWidth', 2);
    hold on
    xlim([-20 100])
    plot([0 0], [-0.1 0.1], '--k', 'LineWidth', 1.5)
    set(gcf, 'Color', [1 1 1])
    ylabel('AEP amplitude (z-score)')
    xlabel('Time (ms)')
    title('AEPs from each rat averaged across all trials')
end
lfp = lfp./length(ratID);
lfp_z = lfp_z./length(ratID);

%% 4. Mean power spectral density in response to 22 kHz single call

% %%%%%%%%%%%%%%%%%%%%%%%% -----  dmPFC ----- %%%%%%%%%%%%%%%%%%%%%%%%


 fileLocation = 'E:\LFPs\AEP2019\Rat';
 
 % Rat IDs for 22kSingleCall 
 
 ratID        = {'26'; '28'; '29'; '30'; '32'; '37'; '39'; '45'; '47'; '48';... 
                 '49'; '50'; '51'; '52'; '53'; '54';'55'; '56'; '59'; '61';...                  
                 '62'; '63'; '64'; '66';'67'; '71'; '73'; '74'; '75'; '76'};
            
                  
 exptID       = '22kSingleCall';
 fileID       = 'SpectralPower';
 brainArea    = 'dmPFC';
 
 % Selected channels for 22kSingleCall

 channelID    = {'31'; '31'; '31'; '31'; '31'; '30'; '30'; '30'; '30'; '29'; '30';...
                 '30'; '30'; '30'; '30'; '30'; '30'; '30'; '30'; '30'; '31';... 
                 '30'; '30'; '30'; '30'; '30'; '30'; '30'; '30'; '30'};
             
%  channelID    = {'32'; '32'; '32'; '32'; '32'; '31'; '31'; '31'; '31'; '30'; '31';...
%                  '31'; '31'; '31'; '31'; '31'; '31'; '31'; '31'; '31'; '32';... 
%                  '31'; '31'; '31'; '31'; '31'; '31'; '31'; '31'; '31'; '31'};
             
%  channelID    = {'31'; '31'; '31'; '31'; '31'; '30'; '30'; '30'; '30'; '29'; '30';...
%                  '30'; '30'; '30'; '30'; '30'; '30'; '30'; '30'; '30'; '31';... 
%                  '30'; '30'; '30'; '30'; '30'; '30'; '30'; '30'; '30'; '30'};

 
 
for i = 1:length(ratID)
    load(strcat(fileLocation, ratID{i, 1}, '\', exptID, '\', fileID, brainArea, '_',  channelID{i, 1}, '.mat'));
end

tf2all = zeros(90, 2001);


% Determining the number of rows and columns in the figure
figure;
nCol = 5;
Q = floor(length(ratID)/nCol); 
R = mod(length(ratID), nCol);

if R == 0
    nRow = Q;
else
    nRow = Q + 1;
end

% Baseline for normalization
baselinetime = [-400 -100];
baselineidx = dsearchn(tx_26', baselinetime');

for i = 1:length(ratID)
    
    tf2 = evalin('base', strcat('tf2_', ratID{i}));
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

    tf2all = tf2all + dbconverted;
%     tf2all = tf2all + pctchange;
%     tf2all = tf2all + baselinediv;
%     tf2all = tf2all + baselineZ;
    subplot(nRow, nCol, i)
    imagesc(tx_54, frex_54, dbconverted); colorbar;
    hold on
    plot([0 0], [1 150], '--k', 'LineWidth', 1);
    set(gca, 'YDir', 'normal', 'CLim', [-1.5 1.5]);
    set(gcf, 'Color', [1 1 1]);
    colormap jet
    ylabel('Frequency (Hz)')
    xlabel('Time (ms)')
    xlim([-400 1400])
    title(strcat('Rat', ratID{i}))
    
end

tf2all = tf2all./length(ratID);
figure;
contourf(tx_54, frex_54, tf2all, 180, 'LineColor', 'none');
hold on
plot([0 0], [1 150], '--k', 'LineWidth', 1.5);
set(gca, 'CLim', [-1.5 1.5])
colormap jet
xlabel('Time(ms)'); ylabel('Frequency (Hz)');
xlim([-400 1400])
title('Power spectral density averaged across all rats and trials')


%% 5. Mean power spectral density in response to 10 kHz pips

% %%%%%%%%%%%%%%%%%%%%%%%% -----  dmPFC ----- %%%%%%%%%%%%%%%%%%%%%%%%


 fileLocation = 'E:\LFPs\AEP2019\Rat';
 
 
 % Rat IDs for 10kpips    
 
  ratID        = {'26'; '28'; '29'; '30'; '32'; '39'; '45'; '47'; '48';... 
                 '49'; '51'; '52'; '53'; '54';'55'; '56'; '59'; '61'; '62';...                  
                 '63'; '64';'71'; '73'; '74'; '75'; '76'; '77'};      
             
                  
 exptID       = '10kpips';
 fileID       = 'SpectralPower';
 brainArea    = 'dmPFC';
 

 
%  Selected channels for 10kpips

 channelID    = {'31'; '31'; '31'; '31'; '31'; '30'; '30'; '30'; '29'; '30';...
                 '30'; '30'; '30'; '30'; '30'; '30'; '30'; '30'; '31';... 
                 '30'; '30'; '30'; '30'; '30'; '30'; '30'};
             
%  channelID    = {'32'; '32'; '32'; '32'; '32'; '31'; '31'; '31'; '30'; '31';...
%                  '31'; '31'; '31'; '31'; '31'; '31'; '31'; '31'; '32';... 
%                  '31'; '31'; '31'; '31'; '31'; '31'; '31'};
% 
% channelID    = {'31'; '31'; '31'; '31'; '31'; '30'; '30'; '30'; '29'; '30';...
%                  '30'; '30'; '30'; '30'; '30'; '30'; '30'; '30'; '31';... 
%                  '30'; '30'; '30'; '30'; '30'; '30'; '30'}; 

for i = 1:length(ratID)
    load(strcat(fileLocation, ratID{i, 1}, '\', exptID, '\', fileID, brainArea, '_',  channelID{i, 1}, '.mat'));
end

tf2all = zeros(90, 2001);


% Determining the number of rows and columns in the figure
figure;
nCol = 5;
Q = floor(length(ratID)/nCol); 
R = mod(length(ratID), nCol);

if R == 0
    nRow = Q;
else
    nRow = Q + 1;
end

% Baseline for normalization
baselinetime = [-400 -100];
baselineidx = dsearchn(tx_54', baselinetime');

for i = 1:length(ratID)
    
    tf2 = evalin('base', strcat('tf2_', ratID{i}));
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

    tf2all = tf2all + dbconverted;
%     tf2all = tf2all + pctchange;
%     tf2all = tf2all + baselinediv;
%     tf2all = tf2all + baselineZ;
    subplot(nRow, nCol, i)
    imagesc(tx_54, frex_54, dbconverted); colorbar;
    hold on
    plot([0 0], [1 150], '--k', 'LineWidth', 1);
    set(gca, 'YDir', 'normal', 'CLim', [-1.5 1.5]);
    set(gcf, 'Color', [1 1 1]);
    colormap jet
    ylabel('Frequency (Hz)')
    xlabel('Time (ms)')
    xlim([-400 1400])
    title(strcat('Rat', ratID{i}))
    
end

tf2all = tf2all./length(ratID);
figure;
contourf(tx_54, frex_54, tf2all, 180, 'LineColor', 'none');
hold on
plot([0 0], [1 150], '--k', 'LineWidth', 1.5);
set(gca, 'CLim', [-1.5 1.5])
colormap jet
xlabel('Time(ms)'); ylabel('Frequency (Hz)');
xlim([-400 400])
title('Power spectral density averaged across all rats and trials')
   


%% end of code
