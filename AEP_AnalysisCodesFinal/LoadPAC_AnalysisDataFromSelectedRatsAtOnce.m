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

 ratID        = {'54'; '57'; '58'; '59'; '62'; ...
                 '64'; '67'; '75'; '81'; '82'}; %             
                  
 exptID       = '22kSingleCall';
% exptID       = '10kpips';
fileID       = 'PACTort';
% brainArea    = 'dmPFC';
brainArea    = 'BLA';

spacer = '_';
 
 
% Selected channels for 22kSingleCall
               
  channelID    = { 'CSC28'; 'CSC28'; 'CSC26'; 'CSC28'; 'CSC28'; ...
                   'CSC28'; 'CSC28'; 'CSC29'; 'CSC01'; 'CSC29'};
 
 

for i = 1:length(ratID)
     
    load(strcat(fileLocation, ratID{i}, '\', exptID, '\', fileID, spacer, brainArea, spacer,...
     ratID{i}(end-1:end), spacer,  channelID{i, 1}, '.mat'));

end

% clear fileLocation ratID exptID fileID brainArea channelID 

%%

Comodulogram_baseline = zeros(length(ratID), length(PhaseFreqVector_CSC28__54), length(AmpFreqVector_CSC28__54));
Comodulogram_stim = zeros(length(ratID), length(PhaseFreqVector_CSC28__54), length(AmpFreqVector_CSC28__54));
MeanAmp_baseline = zeros(10, 1, 36);
MeanAmp_stim = zeros(10, 1, 36);


for i = 1:length(ratID)
    Comodulogram_baseline(i, :, :) = evalin('base', strcat('Comodulogram_baseline_', channelID{i, 1}, spacer, spacer, ratID{i}(end - 1:end)));
    Comodulogram_stim(i, :, :)     = evalin('base', strcat('Comodulogram_stim_', channelID{i, 1}, spacer, spacer, ratID{i}(end - 1:end)));
%     MeanAmp_baseline(i, :, :)      =  evalin('base', strcat('MeanAmp_baseline_', channelID{i, 1}, spacer, spacer, ratID{i}(end - 1:end)));
%     MeanAmp_stim(i, :, :)          =  evalin('base', strcat('MeanAmp_stim_', channelID{i, 1}, spacer, spacer, ratID{i}(end - 1:end)));
%     MI_baseline(i, :)              = evalin('base', strcat('MI_baseline_', channelID{i, 1}, spacer, spacer, ratID{i}(end - 1:end)));
%     MI_stim(i, :)                  = evalin('base', strcat('MI_stim_', channelID{i, 1}, spacer, spacer, ratID{i}(end - 1:end)));
end

Comodulogram_baseline_median = squeeze(median(Comodulogram_baseline, 1));
Comodulogram_stim_median     = squeeze(median(Comodulogram_stim, 1));
% MeanAmp_baseline_median      = squeeze(median(MeanAmp_baseline, 1));
% MeanAmp_stim_median          = squeeze(median(MeanAmp_stim, 1));
% MI_baseline_median           = median(MI_baseline, 1);
% MI_stim_median               = median(MI_stim, 1);

PhaseFreq_BandWidth = 2;
AmpFreq_BandWidth = 20;

%% Plotting
figure;
pcolor(PhaseFreqVector_CSC28__54 + PhaseFreq_BandWidth/2, AmpFreqVector_CSC28__54 + AmpFreq_BandWidth/2, Comodulogram_baseline_median')
shading interp
set(gca, 'FontSize', 14)
ylabel('Amplitude Frequency (Hz)')
xlabel('Phase Frequency (Hz)')
colorbar; colormap jet

figure;
pcolor(PhaseFreqVector_CSC28__54 + PhaseFreq_BandWidth/2, AmpFreqVector_CSC28__54 + AmpFreq_BandWidth/2, Comodulogram_stim_median')
shading interp
set(gca, 'FontSize', 14)
ylabel('Amplitude Frequency (Hz)')
xlabel('Phase Frequency (Hz)')
colorbar; colormap jet

figure;
pcolor(PhaseFreqVector_CSC28__54 + PhaseFreq_BandWidth/2, AmpFreqVector_CSC28__54 + AmpFreq_BandWidth/2, (Comodulogram_stim_median - Comodulogram_baseline_median)')
shading interp
set(gca, 'FontSize', 14)
ylabel('Amplitude Frequency (Hz)')
xlabel('Phase Frequency (Hz)')
colorbar; colormap jet

% figure; clf; 
% bar(10:20/2:720, [ MeanAmp_baseline_median', MeanAmp_baseline_median' ]/sum(MeanAmp_baseline_median'), 'k')
% xlim([0 720])
% set(gca, 'XTick', 0:360:720)
% xlabel('Phase (Deg)')
% ylabel('Amplitude')
% title([ 'MI = ' num2str(MI_baseline_median) ])
% 
% figure; clf; 
% bar(10:20/2:720, [ MeanAmp_stim_median', MeanAmp_stim_median' ]/sum(MeanAmp_stim_median'), 'k')
% xlim([0 720])
% set(gca, 'XTick', 0:360:720)
% xlabel('Phase (Deg)')
% ylabel('Amplitude')
% title([ 'MI = ' num2str(MI_baseline_median) ])


%% Calculating Modulation Index for specific band pairs
% lfp_baseline and lfp_stim need to be defined

Pf1 = 7; % phase frequency 
Pf2 = 12;
Af1 = 80;
Af2 = 100;

[ MI_baseline, MeanAmp_baseline ] = ModIndex_v1(lfp_baseline, Fs, Pf1, Pf2, Af1, Af2, position);
[ MI_stim, MeanAmp_stim ] = ModIndex_v1(lfp_stim, Fs, Pf1, Pf2, Af1, Af2, position);

figure; clf; 
bar(10:20/2:720, [ MeanAmp_baseline, MeanAmp_baseline ]/sum(MeanAmp_baseline), 'k')
xlim([0 720])
set(gca, 'XTick', 0:360:720)
xlabel('Phase (Deg)')
ylabel('Amplitude')
title([ 'MI = ' num2str(MI_baseline) ])

figure; clf; 
bar(10:20/2:720, [ MeanAmp_stim, MeanAmp_stim ]/sum(MeanAmp_stim), 'r')
xlim([0 720])
set(gca, 'XTick', 0:360:720)
xlabel('Phase (Deg)')
ylabel('Amplitude')
title([ 'MI = ' num2str(MI_stim) ])

%% end of code
