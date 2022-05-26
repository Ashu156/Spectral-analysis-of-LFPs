tic;

close all;
clear all;
clc;
%%
%------------ Calcium signal processing from the raw traces file-----------

%% Loading the raw data into MATLAB
[file, path]=uigetfile;
rawDataTraces = readtable(strcat(path,file));
dataTraces = table2cell(rawDataTraces);
headerTraces = (dataTraces(1,2:end));

%% Changing cell IDs into numeric values

nodes = rawDataTraces.Properties.VariableNames(1,2:121);
chrTraces = char(nodes);

for i = 1:size(rawDataTraces,2)-1
    chr2Traces = strsplit(chrTraces(i,:),'C');
    chr3Traces(i) = str2num(chr2Traces{1,2});
end

cellIds = chr3Traces'; % Original cell Ids
cellIds = cellIds +1; % Adding 1 to each cell ID so that cell 'C000' becomes cell 1

%% Extracting signa values from all the cells

valuesTraces = dataTraces(2:end,:); % all the signal values in a table format

% Converting table into a matrix  

nrows = size(valuesTraces,1); % number of rows in the matrix
nCols = size(valuesTraces,2); % number of columns in the matrix
values2Traces = zeros(nrows, nCols); % initiating a zeros matrix 

for i = 1:nrows
    for j  = 1:nCols
        values2Traces(i,j) = str2double(valuesTraces(i,j));
    end
end

timeTraces = values2Traces(:,1); % time vector
SignalTraces = values2Traces(:,2:end); % data from all cells

% Choosing data from those cells which have an accepted tag associated with them 

selectedSignalTraces = []; % data from selected cells
selectedCellids = {}; % selected cell IDs

for k = 1:length(headerTraces)
    token = 'accepted';
    if regexp(token, char(headerTraces(k)),'ignorecase')
        selectedCellids{k} = nodes{k};
        selectedCellids = selectedCellids(~cellfun('isempty', selectedCellids));
        selectedSignalTraces = [selectedSignalTraces SignalTraces(:,k)];
        
    end
end


% Plotting data from all cells
% figure(1);
% surf(timeTraces, [1:numel(cellIds)], SignalTraces'); shading interp; view([118 80]);
% set(gcf,'Color',[1 1 1]); grid off;
% xlabel('Time (s)','FontSize',12); ylabel('Cell ID','FontSize',12); zlabel('Signal Intensity','FontSize',12);
% title('Signal from all cells','FontSize',15)

% Plotting data from accepted cells
figure(2);
surf(timeTraces, 1:numel(selectedCellids), selectedSignalTraces'); shading interp; view([118 80]);
set(gcf,'Color',[1 1 1]); grid off;
xlabel('Time (s)','FontSize',12); ylabel('Cell ID','FontSize',12); zlabel('Signal Intensity','FontSize',12);
title('Signal from accepted cells','FontSize',15);
colors = 'wk';
[cmap] = buildcmap(colors);
colormap(cmap);


%%
%-------------------- Behavioural data processing ------------------------

%% Loading the behavioural events file in  MATLAB

[file, path]=uigetfile;
rawDataEvents = readtable(strcat(path,file));
rawDataEvents = rawDataEvents(16:45,1:3); % limiting to dataEvents of interest, i.e., excluding the headers here
dataEvents = table2cell(rawDataEvents); % converting the table information into a cell array
chrEventTime = char(dataEvents(:,1)); % time information as a string array
eventTimes = str2num(chrEventTime); % converting string into numbers 
state = strcat(dataEvents(:,2),{''},dataEvents(:,3)); % creating all combinations of behavioural states
allStates = unique(state);          % basis vector of behavioural states

%% 
nrowsFinalMatrix = length(state);  % number of rows in the final matrix
nColsFinalMatrix = length(allStates); % number of columns in the final matrix

finalMatrix = zeros(nrowsFinalMatrix, nColsFinalMatrix); % initiating the final matrix as a zeros matrix

%% Color codes for different behavioural states:
% Green: grooming
% Yellow: LED ON
% cyan: exploration in cage
% Blue: halt
% Magenta: walking

colors = {'g', 'y', 'c', 'b', 'm'}; 


%% Assigning START and STOP eventTimes to respective columns of the final matrix for all the distinct behvaioural states

for i = 1:length(state)
    token = state{i}; % gives the first behavioural state found
    
    switch (token) % for details see switch documentation
        
        case char(allStates(1)) % state: START grooming (1st column of final matrix)
            
            if regexp(token, char(allStates(1)), 'ignorecase')             
                chrEventTimes = dataEvents{i,1};
                finalMatrix(i,1) = str2num(chrEventTimes);
                groomingStart(1,i) = str2num(chrEventTimes);
                [~, idx] = find(groomingStart);
                groomingStart = groomingStart(idx); %- eventTimes(1);
            end
        
        case char(allStates(2)) % state: STOP grooming (2nd column of final matrix)
            
        if regexp(token, char(allStates(2)), 'ignorecase')
            chrEventTimes = dataEvents{i,1};
            finalMatrix(i,2) = str2num(chrEventTimes);
            groomingStop(1,i) = str2num(chrEventTimes);
            [~, idx] = find(groomingStop);
            groomingStop = groomingStop(idx);% - eventTimes(1);
        end
        
        for i = 1:length(groomingStart)
            patch([groomingStart(i) groomingStop(i) groomingStop(i) groomingStart(i)],[0 0 (numel(selectedCellids) + 5) (numel(selectedCellids) + 5)],colors{1},'FaceAlpha',0.5, 'EdgeColor', 'none');
        end
        
        case char(allStates(3)) % state: LED ON (3rd column of final matrix)
            
            if regexp(token, char(allStates(3)),'ignorecase')
                chrEventTimes = dataEvents{i,1};
                finalMatrix(i,3) = str2num(chrEventTimes);
                LEDon(1,i) = str2num(chrEventTimes);
                [~, idx] = find(LEDon);
                LEDon = LEDon(idx);% - eventTimes(1);
            end 
            
           
              
            case char(allStates(4)) % state: LED OFF (3rd column of final matrix)
            
            if regexp(token, char(allStates(4)),'ignorecase')
                chrEventTimes = dataEvents{i,1};
                finalMatrix(i,4) = str2num(chrEventTimes);
                LEDoff(1,i) = str2num(chrEventTimes);
                [~, idx] = find(LEDoff);
                LEDoff = LEDoff(idx);% - eventTimes(1);
            end 
            
%              for i = 1:length(LEDon)
%             patch([LEDon(i) LEDoff(i) LEDoff(i) LEDon(i)],[0 0 (numel(selectedCellids) + 5) (numel(selectedCellids) + 5)],colors{1},'FaceAlpha',0.5, 'EdgeColor', 'none');
%              end
        
            case char(allStates(5)) % state: START exploration in cage (1st column of final matrix)
                
                if regexp(token, char(allStates(5)), 'ignorecase')
                    chrEventTimes = dataEvents{i,1};
                    finalMatrix(i,5) = str2num(chrEventTimes);
                    explorationStart(1,i) = str2num(chrEventTimes);
                    [~, idx] = find(explorationStart);
                    explorationStart = explorationStart(idx);% - eventTimes(1);
                end
                
            
            
        
        case char(allStates(6)) % state: STOP exploration in cage (1st column of final matrix)
            
            if regexp(token, char(allStates(6)),'ignorecase')
                chrEventTimes = dataEvents{i,1};
                finalMatrix(i,6) = str2num(chrEventTimes);
                explorationStop(1,i) = str2num(chrEventTimes);
                [~, idx] = find(explorationStop);
                explorationStop = explorationStop(idx);% - eventTimes(1);
            end
            
        
        for i = 1:length(explorationStart)
        patch([explorationStart(i) explorationStop(i) explorationStop(i) explorationStart(i)],[0 0 (numel(selectedCellids) + 5) (numel(selectedCellids) + 5)],colors{3},'FaceAlpha',0.5, 'EdgeColor', 'none');
        end
        
        case char(allStates(7)) % state: START halt/stop (1st column of final matrix)
            
        if regexp(token, char(allStates(7)), 'ignorecase')           
            chrEventTimes = dataEvents{i,1};
             finalMatrix(i,7) = str2num(chrEventTimes);
             haltStart(1,i) = str2num(chrEventTimes);
            [~, idx] = find(haltStart);
            haltStart = haltStart(idx);% - eventTimes(1);
        end
        
        case char(allStates(8)) % state: STOP halt/stop (1st column of final matrix)
            
        if regexp(token, char(allStates(8)), 'ignorecase')            
            chrEventTimes = dataEvents{i,1};
            finalMatrix(i,8) = str2num(chrEventTimes);
            haltStop(1,i) = str2num(chrEventTimes);
            [~, idx] = find(haltStop);
            haltStop = haltStop(idx);% - eventTimes(1);
        end
        
        for i = 1:length(haltStart)
        patch([haltStart(i) haltStop(i) haltStop(i) haltStart(i)],[0 0 (numel(selectedCellids) + 5) (numel(selectedCellids) + 5)],colors{4},'FaceAlpha',0.5, 'EdgeColor', 'none');
        end
        
         case char(allStates(9)) % state: START walking (1st column of final matrix)
            
            if regexp(token, char(allStates(9)), 'ignorecase')          
                chrEventTimes = dataEvents{i,1};
                finalMatrix(i,9) = str2num(chrEventTimes);
                walkingStart(1,i) = str2num(chrEventTimes);
                [~, idx] = find(walkingStart);
                walkingStart = walkingStart(idx);% - eventTimes(1);
            end
        
         case char(allStates(10)) % state: STOP walking (1st column of final matrix)
             
        if regexp(token, char(allStates(10)), 'ignorecase')           
            chrEventTimes = dataEvents{i,1};
            finalMatrix(i,10) = str2num(chrEventTimes);
            walkingStop(1,i) = str2num(chrEventTimes);
            [~, idx] = find(walkingStop);
            walkingStop = walkingStop(idx);% - eventTimes(1);
        end
        
        for i = 1:length(walkingStart)
            patch([walkingStart(i) walkingStop(i) walkingStop(i) walkingStart(i)],[0 0 (numel(selectedCellids) + 5) (numel(selectedCellids) + 5)],colors{5},'FaceAlpha',0.9, 'EdgeColor', 'none');
        end
    end
end

toc;

%% end of code