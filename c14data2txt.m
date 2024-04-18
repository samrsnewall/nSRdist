%% ReadMe
% This script uses my corechoices_MSPF.xlsx spreadsheet and the Mulitza et
% al., 2022 World Atlas data to write all MSPF data from that atlas (that
% fits the criteria I use) into .txt files.

% BEFORE RUNNING, DOUBLE CHECK WHERE THE DATA WILL BE STORED BY LOOKING
% INTO THE FUNCTION netCDF2txt.m!!!!!!!!!

%% Add Functions Folder to Path
addpath('Functions')

%% Load Metadata of MSPF cores
%Check which cores have MSPF (monospecific planktonic foram) dates
data = readtable("../CoreSpreadsheets/COPYcorechoices_MSPF.xlsx"); %read all metadata
dataMSPF = data(data.MSPF == 1,:);

%% Index Good Cores
%Create index for all cores
numAllCores = length(dataMSPF.CoreName);
allcores = 1:numAllCores;
%Create index of cores that have many reversals (determined by manual
%inspection)
reversalDenseCores = ["GeoB1711-4", "H214", "SO75_3_26KL", "MD95-2042"];
badLog = contains(string(dataMSPF.CoreName),reversalDenseCores);
goodLog = badLog == 0;
badIndexes = allcores(badLog);
%Hence create index to use only cores that passed manual inspection
goodIndexes = allcores(~badLog);

%% Choose cores to analyse
%Decide which subset of cores to look at
subsetchooser = zeros(size(goodLog));
subsetchooser(1:5) = 1;
subsetchooserLog = logical(subsetchooser);
chosenCoresLog = goodLog;
cores = table2array(dataMSPF(chosenCoresLog, "CoreName")); %take list of MSPF corenames
lats = table2array(dataMSPF(chosenCoresLog, "LatitudeDec"));
longs = table2array(dataMSPF(chosenCoresLog,"LongitudeDec"));
depths = table2array(dataMSPF(chosenCoresLog, "WaterDepthM"));
LabIDs = table2array(dataMSPF(chosenCoresLog, "LabIDs")); %take list of LabIDs relating to MSPF dates of each core
incDepths = table2array(dataMSPF(chosenCoresLog, "IncludeDepths")); % take list of depths (useful if no labels)
excLabIDs = table2array(dataMSPF(chosenCoresLog, "excludeLabIDs")); %take list of manually removed dates for each core
excDepths = table2array(dataMSPF(chosenCoresLog, "excludeDepth")); %take list of manually removed dates for each core (useful if no labels)
numCores = sum(chosenCoresLog);

%% Output all cores filtered data to txt files
for i = 1:numCores
    netCDF2txt(cores{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i});
end

