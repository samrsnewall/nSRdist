%% ReadMe
% This script uses my corechoices_MSPF.xlsx spreadsheet and the Mulitza et
% al., 2022 World Atlas data to write all MSPF data from that atlas (that
% fits the criteria I use) into .txt files.

% BEFORE RUNNING, DOUBLE CHECK WHERE THE DATA WILL BE STORED!

%% Add Functions Folder to Path
addpath('../Functions')

S.c14AgeLim         = [0 50]; 
S.useLin            = true;
S.usePF             = true;
S.minNumberOfAges   = 4;
S.sheet             = "../DataSheets/COPYcore40MetadataAndLin2014_2.xlsx";
S.folderName = "c14TrainingData_ReversalsRemoved_Feb28";
S.removeLargeGaps = false;
S.WApath       = "/Applications/PaleoDataView/WA_Foraminiferal_Isotopes_2022";
S.sandboxPath  = "/Volumes/ExtDrive850X/MATLAB/nSRdist_code";

%% Load Metadata of MSPF cores
%Check which cores have MSPF (monospecific planktonic foram) dates
data = readtable(S.sheet); %read all metadata

if S.useLin == 1 && S.usePF == 1
dataMSPF = data(data.WAuse == 1 | data.Lin2014Keep == 1,:);
end

d = load("../Results/dataT_RLGfalse_R200M20_Feb26.mat");
dataMSPF1 = d.dataT;

%% Index Good Cores
%Create index for all cores
numAllCores = length(dataMSPF.CoreName);
allcores = 1:numAllCores;
%Create index of cores that have many reversals (determined by manual
%inspection)
reversalDenseCores = ["PlaceholderString"];%["GeoB1711-4", "GIK17940-2", "H214", "SO75_3_26KL", "MD95-2042", "RC11-83"];
badLog = contains(string(dataMSPF.CoreName),reversalDenseCores);
restrictionLog =  dataMSPF1.meanSR > 8 & dataMSPF1.depths > 1000 & abs(dataMSPF1.lats) < 40;
goodLog = badLog == 0 & restrictionLog == 1;
badIndexes = allcores(badLog);
%Hence create index to use only cores that passed manual inspection
goodIndexes = allcores(~badLog);

%% Choose cores to analyse
%Decide which subset of cores to look at
subsetchooser = zeros(size(goodLog));
subsetchooser(1:5) = 1;
subsetchooserLog = logical(subsetchooser);
chosenCoresLog = goodLog;
cores = table2array(dataMSPF1(chosenCoresLog, "cores")); %take list of MSPF corenames
lats = table2array(dataMSPF1(chosenCoresLog, "lats"));
longs = table2array(dataMSPF1(chosenCoresLog,"longs"));
depths = table2array(dataMSPF1(chosenCoresLog, "depths"));
[LabIDs, incDepths, excLabIDs, excDepths, dataLoc] = extract3(dataMSPF, chosenCoresLog, S);
numCores = sum(chosenCoresLog);

% %%
% for i = 1:numCores% find(namedLog)'
%   corePlotCal(cores{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, dataLoc(i), S)
% end

%% Output all cores filtered data to txt files
%Decide if I want all cores or only cores with SR > 8 (0 or 1 respectively)
sepBySRg8 = 0;
%Run through all cores to output data to txt file
for i = 1:numCores
    netCDF2txt(cores{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, dataLoc(i),  sepBySRg8, S);
end