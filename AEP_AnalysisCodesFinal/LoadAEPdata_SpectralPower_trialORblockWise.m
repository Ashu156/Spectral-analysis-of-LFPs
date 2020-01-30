%% 0.

clear all; close all; clc;

%% 1. Rat IDs for 22kSingleCall BLA

fileLocation = 'E:\LFPs\AEP2019\Rat';

% Rat IDs for 22kSingleCall

 ratID        = {'54'; '57'; '58';'59'; '62'; ...
                 '64'; '67'; '75';'81'; '82'}; %             
                  
 exptID       = '22kSingleCall';
 fileID       = 'SpectralPower';
 brainArea    = 'BLATrialWise';
 nCycles      = '40cycles'; 
 
% Selected channels for 22kSingleCall
               
  channelID    = { '28'; '28'; '26'; '28'; '28'; ...
                  '28'; '28'; '28'; '01'; '29'};
 
 

for i = 1:length(ratID)   % 
    load(strcat(fileLocation, ratID{i, 1}, '\', exptID, '\', fileID, brainArea, '_',  channelID{i, 1}, '_', nCycles, '_1000ms', '.mat'));
% load(strcat(fileLocation, ratID{i, 1}, '\', exptID, '\', fileID, brainArea, '_',  channelID{i, 1},  '_', nCycles, '_100ms', '.mat'));
end

%%
baseline = [-1000 0];
baselinetimeidx = dsearchn(tx', baseline');

for i = 1:length(ratID)
    
    tf = evalin('base', strcat('tf_', ratID{i}));
    baseline_powerZ = tf(:, baselinetimeidx(1):baselinetimeidx(2), :);
    baselineZ{i, 1} = (tf - repmat(mean(baseline_powerZ, 2), 1, size(tf, 2))) ./ repmat(std(baseline_powerZ, [], 2), 1, size(tf, 2));
%     tf_10(:, :, i)  = mean(tf(:, :, 1:10), 3);
%     tf_20(:, :, i)  = mean(tf(:, :, 11:20), 3);
%     tf_30(:, :, i)  = mean(tf(:, :, 21:30), 3);
%     tf_40(:, :, i)  = mean(tf(:, :, 31:40), 3);
%     tf_50(:, :, i)  = mean(tf(:, :, 41:50), 3);
%     tf_60(:, :, i)  = mean(tf(:, :, 51:60), 3);
%     tf_70(:, :, i)  = mean(tf(:, :, 61:70), 3); 
%     tf_80(:, :, i)  = mean(tf(:, :, 71:80), 3); 
%     tf_90(:, :, i)  = mean(tf(:, :, 81:90), 3);
%     tf_100(:, :, i) = mean(tf(:, :, 91:100), 3);
    
end

for i = 1:length(ratID)
    
%     fprintf('Running for Rat #%d\n', i);
    for k = 1:100
        
        fprintf(['Running Trial #%d\n for Rat %d\n'], k, i)
        trial = baselineZ{i, 1}(:, :, k);
        
        % New variable names
    newname = strcat('trial', num2str(k),'_',  ratID{i});
    
    % Assigning the old values to the renamed variables
    str = [newname, ' = ', 'trial', ';'];
    
    
        evalin('base',str)
   
         save(strcat(fileLocation, ratID{i}, '\', exptID,'\', 'CSC', channelID{i, 1}, '_', 'trial', num2str(k), '_', 'zScoredPower_BLA', '_', nCycles, '_1000ms'),strcat('trial',  num2str(k), '_', num2str(ratID{i})))

    end
    
end