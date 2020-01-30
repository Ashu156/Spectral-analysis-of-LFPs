%%
clear all; close all; clc;

%%
trialId = 3;

load(strcat('E:\LFPs\AEP2019\trial', num2str(trialId),'_AllRats_zScoredPower_BLA_40cycles_1000ms.mat'));

tx = [-1500:1:2500];
baseline = [-1000 0];
stim1     = [1 1300];
stim2     = [25 1300];


theta      = [4 12];
slowGamma   = [30 70];
fastGamma  = [70 120];

baselinetimeidx = dsearchn(tx', baseline');
stim1timeidx    = dsearchn(tx', stim1');
stim2timeidx    = dsearchn(tx', stim2');

%% end of code

for i = 1:10
    trial = evalin('base', strcat('trial', num2str(trialId)));
    
     baseline_data = trial{1,i}(:,baselinetimeidx(1):baselinetimeidx(2));
     baseline_meandata = mean(baseline_data, 2);
     baseline_theta(i, 1) = sum(baseline_meandata(theta(1):theta(2)));
     baseline_slowGamma(i, 1) = sum(baseline_meandata(slowGamma(1):slowGamma(2)));
     baseline_fastGamma(i, 1) = sum(baseline_meandata(fastGamma(1):fastGamma(2)));
     
     stim1_data = trial{1,i}(:,stim1timeidx(1):stim1timeidx(2));
     stim1_meandata = mean(stim1_data, 2);
     stim1_theta(i, 1) = sum(stim1_meandata(theta(1):theta(2)));
     stim1_slowGamma(i, 1) = sum(stim1_meandata(slowGamma(1):slowGamma(2)));
     stim1_fastGamma(i, 1) = sum(stim1_meandata(fastGamma(1):fastGamma(2)));
     
     stim2_data = trial{1,i}(:,stim2timeidx(1):stim2timeidx(2));
     stim2_meandata = mean(stim2_data, 2);
     stim2_theta(i, 1) = sum(stim2_meandata(theta(1):theta(2)));
     stim2_slowGamma(i, 1) = sum(stim2_meandata(slowGamma(1):slowGamma(2)));
     stim2_fastGamma(i, 1) = sum(stim2_meandata(fastGamma(1):fastGamma(2)));
 end