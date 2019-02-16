% This piece of code processes the data recorded from Neuralynx's data acquisition system 
% and saves your data of interest in a .mat file for further anlysis

% loading .nse spike files into MATLAB
clear all; close all; clc;

%%-----------AEP------------------------------
% Get directory contents
direct='F:\LFPs\AEP2019\Rat29\10kpips\'; % source directory for files
dirData=dir(direct);
dirData = dirData(~[dirData.isdir]);  % Use only the file data
fileNames = {dirData.name};           % Get file names

nse = '.nse';
nev = '.nev';
ncs = '.ncs';
nvt = '.nvt';
fileData_nse = regexp(fileNames,nse,'tokens');       %# Find tokens
fileData_nev = regexp(fileNames,nev,'tokens');
fileData_ncs = regexp(fileNames,ncs,'tokens');
fileData_nvt=regexp(fileNames,nvt,'tokens');
index_nse = ~cellfun('isempty',fileData_nse);        %# Find index of matches
index_nev = ~cellfun('isempty',fileData_nev);
index_ncs = ~cellfun('isempty',fileData_ncs);
index_nvt=~cellfun('isempty',fileData_nvt);
fileNames_nse = fileNames(index_nse);                %# Remove non-matching file names
fileNames_nev = fileNames(index_nev);
fileNames_ncs = fileNames(index_ncs);
fileNames_nvt=fileNames(index_nvt);

% % %loading *.nse files
for i=1:length(fileNames_nse)
      file_nse=char(fileNames_nse(i));
      file_sp=strcat(direct,file_nse);
 
      FieldSelection=[1 0 1 0 1];
      HeaderExtraction=1;
      ExtractMode=1; %Extract All
      [ts_spike,n_cell,dp_spike,header]= Nlx2MatSpike(file_sp, FieldSelection,HeaderExtraction, ExtractMode,[]);
      
      %sorting the cells
      num=max(n_cell);
      for j=1:num
          nn_ts=char(strcat(file_nse(1:5),'_','0',num2str(j),'_','TS'));
          nn_dp=char(strcat(file_nse(1:5),'_','0',num2str(j),'_','DP'));
          assignin('base',nn_ts,ts_spike(n_cell==j));
          assignin('base',nn_dp,dp_spike(:,:,n_cell==j));
      end
      nn_hd=char(strcat(file_nse(1:4),'_','NlxHeader'));
      assignin('base',nn_hd,header);      
end
% % 
%loading *.ncs files
for i=1:length(fileNames_ncs)
      file_ncs=char(fileNames_ncs(i));
      file_csc=strcat(direct,file_ncs);
 
      FieldSelection=[1 0 0 0 1];
      HeaderExtraction=1;
      ExtractMode=1; %Extract All
      [ts_csc,dp_csc,header_csc]= Nlx2MatCSC(file_csc, FieldSelection,HeaderExtraction, ExtractMode);
      
      %sorting the cells
      nn_csc_ts=char(strcat(file_ncs(1:5),'_','TS'));
      nn_csc_dp=char(strcat(file_ncs(1:5),'_','DP'));
      nn_hd_csc=char(strcat(file_ncs(1:5),'_','NlxHeader'));
      assignin('base',nn_csc_ts,ts_csc);
      assignin('base',nn_csc_dp,dp_csc);
      assignin('base',nn_hd_csc,header_csc);
end

% % % loading *.nev files
for i=1:length(fileNames_nev)
      file_nev=char(fileNames_nev(i));
      file_ev=strcat(direct,file_nev);
 
      FieldSelection=[1 0 1 0 1];
      HeaderExtraction=0;
      ExtractMode=1; %Extract All
      [ts_events,ttl_ev,string_ev]= Nlx2MatEV(file_ev, FieldSelection,HeaderExtraction, ExtractMode,[]);
      
      %sorting the events
      nn_ev=char(strcat(file_nev(1:6)));
      nn_ev_ttl=char(strcat(file_nev(1:6),'_','ttl'));
      nn_ev_string=char(strcat(file_nev(1:6),'_','string'));
      assignin('base',nn_ev,ts_events);
      assignin('base',nn_ev_ttl,ttl_ev);
      assignin('base',nn_ev_string,string_ev);
      
      %sorting the events
      % k_1s=find(Events_ttl==-32768);
      % k_p5s=find(Events_ttl==16384);
      % k_p2s=find(Events_ttl==8192);
      % Events_1s=Events(k_1s);
      % Events_p5s=Events(k_p5s);
      % Events_p2s=Events(k_p2s);
      
%       start=find(Events_ttl==4);
      play=find(Events_ttl==4);

%         start=find(Events_ttl==4096*1);
        Events_tone_on=Events(play);
end

%loading *.nvt files
for i=1:length(fileNames_nvt)
      file_nvt=char(fileNames_nvt(i));
      file_coord=strcat(direct,file_nvt);
 
      FieldSelection=[1 1 1 1 1 1];
      HeaderExtraction=1;
      ExtractMode=1; %Extract All
      [ts_coord,extracted_x,extracted_y,extracted_angle,target, point, header_extract]= Nlx2MatVT(file_coord, FieldSelection,HeaderExtraction, ExtractMode,[]);
end
% 
%  clear Extract*; clear Field*; clear Header*; clear dir*; clear file*; clear n*; clear ttl*;clear k_*;
%  clear ts*; clear index*; clear header*; clear dp*; clear string*; clear i; clear j;
% %  clear Events_string; clear Events_ttl;

% plot(extracted_x,extracted_y,'r')
save('F:\LFPs\AEP2019\Rat29\10kpips\matfile.mat')



