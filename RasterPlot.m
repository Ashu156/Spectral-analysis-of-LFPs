close all; 
clear all; 
clc;

%% Loading the data file into MATLAB

[file, path]=uigetfile;        % graphic user interface for selecting the data file
rawData  = readtable(strcat(path,file)); % reading the data from the file into MATLAB
rawData = sortrows(rawData,'Time_s_','ascend'); % sorting values according to time
data = table2cell(rawData); % converting table into cell array

%% Changing cell ids into numeric values

chr = char(data(:,2)); 
for i = 1:size(rawData,1)
    chr2 = strsplit(chr(i,:),'C');
chr3(i) = str2num(chr2{1,2});
end

%%
times = cell2mat(data(:,1)); % defining time values
dt = times(2)-times(1); % sampling period
dataModified (:,1) = cell2mat(data(:,1));
cellIds = chr3'; % Original cell Ids
cellIds = cellIds +1; % Adding 1 to each cell Id so that cell 'C000' becomes cell 1
values = cell2mat(data(:,3)); % Readings from each cell (probably, (delta F/ F) in this case) 
firstCell = min(chr3); % first cell Id
lastCell = max(chr3); % last cell Id
allCells = (firstCell + 1):1:(lastCell + 1); % creating a vector of cell Ids from first to last

%%
nRows = size(allCells,2); % number of rows in the final matrix
nCol = size(times,1);     % number of columns in the final matrix
tf = zeros(nRows, nCol); % initiating the final zero matrix with number of rows equal to number of cells and number of columns equal to number of time points

%% Updating the values of the final matrix

j = 1; % initiating a counter 


for k = 1:size(rawData,1)
tf(cellIds(k),j) = values(k); % updating the values in the final matrix as per the data
j = j+1; % increaing the counetr by unit step
end
%% Plotting the final data
figure('Color',[1 1 1]);surf(times, allCells, tf); colorbar; view([118 80]); shading interp; grid off;
colors = 'wbr';
[cmap] = buildcmap(colors);
colormap(cmap); % colormap of choice
caxis([0 1600]); % limits of the color axis
xlabel('time (s)'); ylabel ('Cell #'); zlabel('Signal intensity') % X-,Y- and Z- axis labels

%% end
