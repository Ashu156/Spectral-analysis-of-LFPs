% This piece of code uses a family of complex Morlet wavelets to decompose
% the signal (Auditory Evoked Potentials, here) in the time-frequency plane


clear all; % clearing the workspace variables
close all; % closing any tab if open from previous execution
clc;       % clearing the command prompt space 

load('F:\LFPs\AEP2019\Rat29\10kpips\matfile.mat'); % loading the .mat file of interest


x_min = -0.5; % decides the minimum time point of consideration (in seconds) of the AEP from event onset 
x_max = 0.5;  % decides the maximum time point of consideration (in seconds) of the AEP from event onset

% CHOOSING THE REQUIRED CHANNEL

ts_csc = CSC31_TS; % Channel of choice
dp_csc = (CSC31_DP); % Data points in channel of choice
Fs = str2num(CSC28_NlxHeader{14}(end-3:end)); % Sampling frequency
Fnyq = round(Fs/2);                           % Nyquist frequency
ADBitVolts = (str2double(CSC31_NlxHeader{15}(13:38)))*10^6; % conversion factor in microVolts
dp_csc = dp_csc*ADBitVolts; % data points in microVolts
dp_csc = dp_csc(:); % linearizing the data points in a single column vector

% dp_csc = detrend(dp_csc,'linear');
ts_events = (Events_tone_on + 000000.00)./10^6; % Selecting events of choice (onset of auditory stimulus in this case)
numTrials = length(ts_events); % number of events
ts_csc = ts_csc./10^6 ; % time-stamps in seconds

% Generate all the timestamps corresponding to all data points
tt = [0:1/Fs:511/Fs]' ;  
tts_csc = [];           % Initializing a new empty matrix to interpolate time-stamps of data points

for i = 1:length(ts_csc)
    tts_csc = [tts_csc  [ts_csc(i) + tt]];
end

tts_csc = tts_csc(:); % Linearizing the interpolated time-stamp series in a single column vector
kmin = round(Fs*(x_min));         % Minimum x-limit  
kmax = round(Fs*(x_max));         % Maximum x-limit  

tx = (x_min:1/Fs:x_max)*1000 ; % time of interest in milliseconds

 

bin_csc_BLA = zeros(length([kmin:kmax]),numTrials); % Initializing a new zero matrix for storing data values of interest where data is arranged as follows: rows represent data values for selected time of interest, whilecolumns represent trial numbers

for i = 1:numTrials
    min_dist(1,i) = min(abs(tts_csc(:)-ts_events(i)));
    k0(1,i) = find(abs(tts_csc(:)-ts_events(i)) == min_dist(1,i)) ;
    bin_csc_BLA(:,i) = dp_csc(k0(1,i)+kmin : k0(1,i)+kmax);
end

bin_csc_BLA = bin_csc_BLA(1:length(tx),:);

%% Plotting the perievent hsitogram
% figure; clf
% subplot(3,1,2:3); hold on
mm = mean(bin_csc_BLA,2);
mm = smooth(mm,0.01*length(bin_csc_BLA),'loess');
   
%% Filtering the signal using a bandpass filter


F_cutL = 2; F_cutH = 60; % Low and High cut off frequencies in Hz
[zz,pp,kk] = ellip(20, 0.2, 80, [F_cutL F_cutH]./Fnyq); % a 20th orderelliptic digital bandpass filter  
[sos,g] = zp2sos(zz,pp,kk);	      % Convert to SOS form
Hd = dfilt.df2tsos(sos,g);        % Create a dfilt object
Data_filt  = filtfilthd(Hd, mm); % Filtered signal
mm = Data_filt;



baseline = [-500 0 500]; % time range of interest in milliseconds

% convert baseline from milliseconds to indices
baseidx = dsearchn(tx',baseline');

baselineData = mean(mm(baseidx(1):baseidx(2)));
normalizedData = mm(:) - baselineData;

%%
% wavelet parameters
min_freq = 1;
max_freq = 150;
num_frex = 100;


% other wavelet parameters
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
% frex = min_freq:0.1:max_freq;
time = -1:1/Fs:1;
half_wave = (length(time)-1)/2;

% FFT parameters
nKern = length(time);
nData = length(normalizedData)*numTrials;
nConv(1:2) = nKern + nData - 1;
nConv(3)   = nKern + length(normalizedData)-1; % ERP is only one trial-length

% initialize output time-frequency data
tf = zeros(length(frex),length(normalizedData));

%% prepare non-phase-locked activity

% compute ERP
erp = normalizedData;

erp_orig = polyfit(tx',erp,1);
erp_detrend = erp - polyval(erp_orig,tx');
erp_dt = polyfit(tx',erp_detrend,1);
figure(1)
plot(tx',erp,'-b',tx', polyval(erp_orig,tx'),'--r')
hold on
plot(t',erp_detrend,'-k',t', polyval(erp_dt,tx'),'--r','linew',2)
hold off
legend('Original data', 'Trend', 'Detrended data', 'Mean of detrended data');

% compute non-phase-locked power by subtracting ERP from each trial
nonphase_EEG = squeeze( bsxfun(@minus,bin_csc_BLA(:,:),erp) );

% for i=1:numTrials
% nonphase_EEG(:,i)=nonphase_EEG(:,i)-mean(nonphase_EEG,2);
% end

figure(2), clf
subplot(311)
plot(tx,erp)
xlabel('Time (ms)'), ylabel('Voltage (\muV)')
% set(gca,'xlim',[-300 1300])

subplot(312)
plot(tx,bin_csc_BLA(:,10))
hold on
plot(tx,squeeze(nonphase_EEG(:,10)),'r')
legend({'total';'non-phase-locked'})
xlabel('Time (ms)'), ylabel('Voltage (\muV)')
% set(gca,'xlim',[-300 1300])

subplot(313)
plot(tx,erp)
hold on
plot(tx,mean(nonphase_EEG,2),'r')
legend({'total ERP';'non-phase-locked ERP'})
xlabel('Time (ms)'), ylabel('Voltage (\muV)')
% set(gca,'xlim',[-300 1300])

%% time-frequency decomposition

% FFT of total data
fft_EEG{1} = fft( reshape(bin_csc_BLA(:,:),1,[]), nConv(1));

% % FFT of ERP (phase-locked data)
% fft_EEG{3} = fft( erp ,nConv(3));

% range_cycles=[2 10];
% num_cycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),50);
% num_cycles = length(range_cycles);

% for cyclei=1:length(num_cycles)
    for fi=1:length(frex)
    
    % create wavelet and get its FFT
    s        = 6 /(2*pi*frex(fi));

    wavelet  = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));
    
    % run convolution for each of total, induced, and evoked
   
%         figure (4); plot(time,wavelet); hold on;
        % need separate FFT 
        waveletX = fft(wavelet,nConv(1));
        waveletX = waveletX./max(waveletX);
        
%         figure (3); plot(frex,abs(waveletX(1:length(frex))));hold on;
        
        % notice that the fft_EEG cell changes on each iteration
        as = ifft(waveletX.*fft_EEG{1},nConv(1));
        as = as(half_wave+1:end-half_wave);
        as = reshape(as,length(normalizedData),numTrials);
        
            % compute power
            temppow = mean(abs(as).^2,2);
        
        
        % db correct power
        tf(fi,:) = 10*log10(temppow);
        tfBaseline(fi,:) = mean(tf (fi,(baseidx(1):baseidx(2))));
        tfEvoked(fi,:) = bsxfun(@minus,tf(fi,:),tfBaseline(fi,:));
%         tfEvoked(fi,:) = tf(fi,:) - repmat(tfBaseline(fi,:),length(tx));
        
     % end loop around total, evoked, induced
end % end frequency loop
% end % end cycle loop


% tfEvoked = tf(:,:) - (repmat(tfBaseline(:,:), length(tx)))';


% color limits
clims = [-2 4; 0 5; -15 15; 0 0.8 ];



figure(2), clf

    
  surf(tx,frex,tfEvoked);view([0 90]); shading interp;
  
set(gca,'clim',clims(1,:),'xlim',[x_min*1e3 x_max*1e3],'xtick',x_min*1e3:100:x_max*1e3,'YDir','normal')
    xlabel('Time (ms)')
    ylabel('Frequency (Hz)')
   



% subplot(235)
% 
% % estimate phase-locked part of the signal as the difference 
% % between the phase-locked and non-phase-locked
% phaselockedTF = squeeze(tf(1,:,:) - tf(2,:,:));
% % the next line is equivalent to the previous line, just FYI
% % phaselockedTF = squeeze(diff(tf([2 1],:,:),[],1));
% 
% contourf(tx,frex,phaselockedTF,40,'linecolor','none')
% set(gca,'clim',clims(1,:),'xlim',[-300 500],'xtick',-200:100:500)
% xlabel('Time (ms)'), ylabel('Frequency (Hz)')
% title('Phase-locked')

% plot ERP on top
% hold on
% plot(tx,erpt,'k','Linewidth',2)

%%
% Calcualting power in different frequency bands during baseline and tone
% ON period

frequencyBand={[25 40]};
% frequencyName={'Theta','Beta','Low Gamma','High Gamma'};

PowerBaseline=cell(1,length(frequencyBand));
PowerEvoked=cell(1,length(frequencyBand));


for i=1:length(frequencyBand)
    freqLow(i) = frequencyBand{i}(1);                                         % finding the lower limits of frequency bands
    freqHigh(i)  = frequencyBand{i}(2);                                       % finding the upper limits of frequency bands
    
    freqLowMinDist(1,i) = min(abs(frex(:)-freqLow(i)));                         % calculating the minimum of the absolute difference between the frequencies and lower limits
    freqLowidx(1,i) = find(abs(frex(:)-freqLow(i)) == freqLowMinDist(1,i)) ;    % finding the indices of lower limits in the frequency vector
    
    freqHighMinDist(1,i) = min(abs(frex(:)-freqHigh(i)));                       % calculating the minimum of the absolute difference between the frequencies and upper limits
    freqHighidx(1,i) = find(abs(frex(:)-freqHigh(i)) == freqHighMinDist(1,i)) ; % finding the indices of upper limits in the frequency vector
    
    Power{1,i} = tfEvoked(freqLowidx(i):freqHighidx(i),baseidx(1):baseidx(3));
    PowerBaseline{1,i} = tfEvoked(freqLowidx(i):freqHighidx(i),baseidx(1):baseidx(2));
    PowerEvoked{1,i} = tfEvoked(freqLowidx(i):freqHighidx(i),baseidx(2):baseidx(3));
    
end
