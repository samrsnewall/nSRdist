%% oneCoreRestrictTest.m
% This script takes one core and calculates the distribution of
% normalised Sed Rates. It does this multiple times, each with a new
% restriction on the minimum allowed age difference

%% Add folder of necessary functions to path
addpath('Functions')

%% Load Metadata of MSPF cores
%Check which cores have MSPF (monospecific planktonic foram) dates
data     = readtable("COPYcorechoices_MSPF_highRes2.xlsx"); %read all metadata
dataMSPF = data(data.MSPF == 1,:);

%% Choose core to look at
%Create index for a specific core
namedLog = contains(string(dataMSPF.CoreName), "SHAK06-5K");

%% Get material data from Excel
% Get info from Excel Spreadsheet
chosenCoresLog = namedLog;
cores       = table2array(dataMSPF(chosenCoresLog, "CoreName")); %take list of MSPF corenames
lats        = table2array(dataMSPF(chosenCoresLog, "LatitudeDec"));
longs       = table2array(dataMSPF(chosenCoresLog,"LongitudeDec"));
depths      = table2array(dataMSPF(chosenCoresLog, "WaterDepthM"));
LabIDs      = table2cell(dataMSPF(chosenCoresLog, "LabIDs")); %take list of LabIDs relating to MSPF dates of each core
incDepths   = table2cell(dataMSPF(chosenCoresLog, "IncludeDepths")); % take list of depths (useful if no labels)
excLabIDs   = table2cell(dataMSPF(chosenCoresLog, "excludeLabIDs")); %take list of manually removed dates for each core
excDepths   = table2cell(dataMSPF(chosenCoresLog, "excludeDepth")); %take list of manually removed dates for each core (useful if no labels)
numCores    = sum(chosenCoresLog);

%% Plot radiocarbon data of chosen core
% Plot all radiocarbon data
iPlot = 1;
corePlot(cores{iPlot}, LabIDs{iPlot}, incDepths{iPlot}, excLabIDs{iPlot}, excDepths{iPlot})

%% Find reversals in the core data
%Get reversal data, etc
i = 1;
[core_invSRvals{i}, core_invSRprobs{i}, meanSR(i), MSI_byage(i), MSI_bydepth(i), sedimentlength(i), num14cpairs(i), corescenarios{i}, newlabels{i}, numreversals(i)] = oneCoreSRpdf(cores{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, 0);

%% Use random sampling approach to get normalised SR data
%Restriction = agediff > 0y
ii = 1;
[nSRcounts{ii}, agediffs{ii}] = oneCoreTMRestrict(cores{i}, corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, 0);
plotSRandResHistograms(nSRcounts(ii), agediffs(ii), num14cpairs, 1, 101, 'k', cores{1})
[~,~,TM0] = TMcalculation(nSRcounts{ii});

%% Use random sampling approach to get normalised SR data with some restriction on acceptable age differences
% Restriction = agediff > 500y
ii = 2;
[nSRcounts{ii}, agediffs{ii}] = oneCoreTMRestrict(cores{i}, corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, 500);
plotSRandResHistograms(nSRcounts(ii), agediffs(ii), num14cpairs, 1, 101, 'k', cores{1})
[~,~,TM500] = TMcalculation(nSRcounts{ii});

% Restriction = agediff > 1000y
ii = 3;
[nSRcounts{ii}, agediffs{ii}] = oneCoreTMRestrict(cores{i}, corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, 1000);
plotSRandResHistograms(nSRcounts(ii), agediffs(ii), num14cpairs, 1, 101, 'k', cores{1})
[~,~,TM1000] = TMcalculation(nSRcounts{ii});

%Restriction = agediff > 1500y
ii = 3;
[nSRcounts{ii}, agediffs{ii}] = oneCoreTMRestrict(cores{i}, corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, 1500);
plotSRandResHistograms(nSRcounts(ii), agediffs(ii), num14cpairs, 1, 101, 'k', cores{1})
[~,~,TM1500] = TMcalculation(nSRcounts{ii});


outputMetadataAndSummaryFigures(1, cores, lats, longs, depths, meanSR, MSI_byage, MSI_bydepth, nSRcounts(1), sedimentlength, num14cpairs)
