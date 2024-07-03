%% Add folder of necessary functions to path
addpath('Functions')

%% Load Metadata of MSPF cores
%Check which cores have MSPF (monospecific planktonic foram) dates
data     = readtable("COPYcorechoices_MSPF_highRes2.xlsx"); %read all metadata
dataMSPF = data(data.MSPF == 1,:);

%% Choose cores to look at
%------ Index Good Cores
%Create index for all cores
numAllCores = length(dataMSPF.CoreName);
allcores    = 1:numAllCores;

%Create logical of cores without those that have been manually removed
%(goodLog)
reversalDenseCores = ["GeoB1711-4", "H214", "SO75_3_26KL", "KNR159-5-36GGC"];
problemCores       = [];
badLog             = contains(string(dataMSPF.CoreName),[reversalDenseCores, problemCores]);
badIndexes         = allcores(badLog);
goodLog            = badLog == 0;
goodIndexes = allcores(~badLog);

%Create some other useful logical, to test all cores, only 1 core (defined by name), or a
%subset of cores (defined by index)
namedLog           = contains(string(dataMSPF.CoreName), "GEOFARKF13");

allLog             = true(length(goodLog), 1);

subsetChooser       = logical(zeros(numAllCores, 1)); %#ok<LOGL>
subsetChooser(1:5)  = 1; 

%% Get material data from Excel
% Get info from Excel Spreadsheet
chosenCoresLog = goodLog;
cores       = table2array(dataMSPF(chosenCoresLog, "CoreName")); %take list of MSPF corenames
lats        = table2array(dataMSPF(chosenCoresLog, "LatitudeDec"));
longs       = table2array(dataMSPF(chosenCoresLog,"LongitudeDec"));
depths      = table2array(dataMSPF(chosenCoresLog, "WaterDepthM"));
LabIDs      = table2cell(dataMSPF(chosenCoresLog, "LabIDs")); %take list of LabIDs relating to MSPF dates of each core
incDepths   = table2cell(dataMSPF(chosenCoresLog, "IncludeDepths")); % take list of depths (useful if no labels)
excLabIDs   = table2cell(dataMSPF(chosenCoresLog, "excludeLabIDs")); %take list of manually removed dates for each core
excDepths   = table2cell(dataMSPF(chosenCoresLog, "excludeDepth")); %take list of manually removed dates for each core (useful if no labels)
numCores    = sum(chosenCoresLog);

%% invSR PDF Approach
% This section runs through a quick, less-computationally-expensive
% estimator of SR to find any reversals, create scenarios, to calculate the
% meanSR and the number of dates used overall

%------ Initialise variables to hold this information
core_invSRvals  = cell(numCores,1);
core_invSRprobs = cell(numCores,1);
meanSR          = nan(numCores,1);
MSI_byage       = nan(numCores,1);
MSI_bydepth     = nan(numCores,1);
sedimentlength  = nan(numCores,1);
num14cpairs     = nan(numCores,1);
transprobs_cores= nan(3,3,numCores);
corescenarios   = cell(numCores,1);
newlabels       = cell(numCores,1);
numreversals    = nan(numCores,1);

%----- Apply invSR with loop
%Calculate SR distribution for each core, as well as meanSR and other
%useful information
for i = 1:numCores
    [core_invSRvals{i}, core_invSRprobs{i}, meanSR(i), MSI_byage(i),...
        MSI_bydepth(i), sedimentlength(i), num14cpairs(i),...
        corescenarios{i}, newlabels{i}, numreversals(i)]...
        = oneCoreSRpdf(cores{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, 0);
end

%% ------ Define Subsets of interest
% high SR and low SR cores (separated by 8cm/kyr following Lin et al., 2014)
lowSRCoresLog   = meanSR<= 8;
highSRCoresLog  = meanSR >8;
allCoresLog     = ~isnan(meanSR);

% find 10 highSR cores with most data
[~,highSRhighResCoresInd]   = maxk(num14cpairs.*highSRCoresLog, 10);
highSRhighResCoresLog       = unfind(highSRhighResCoresInd, numel(cores));

%% invSR Random Sampling Approach
%This section runs a slower, more computationally-expensive estimator of SR
%to construct the SR distributions. The benefits that this method has over
%the previous method is that it can be used to create a transition matrix
%and it can be used to weight by age, and we can create distributions of
%age difference between radiocarbon pairs to test for a resolution effect.

%------ Run through all chosen cores with random sampling approach,
% RESTRICTION ON MINIMUM AGE DIFFERENCE = 0

%Initiate variables
nSRcounts       = cell(numCores, 1); % Holds all the nSR counts (which form histogram that makes nSR pdf)
agediffs        = cell(numCores, 1); % Holds all the age differences for each nSR measurement (resolution pdf)

parfor i =  1:numCores
    [nSRcounts{i}, agediffs{i}] = oneCoreTMRestrict(cores{i}, corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, 0);
end

%------ Run through all chosen cores with random sampling approach,
% RESTRICTION ON MINIMUM AGE DIFFERENCE = 500

%Initiate variables
nSRcounts500   = cell(numCores, 1); % Holds all the nSR counts (which form histogram that makes nSR pdf)
agediffs500    = cell(numCores, 1); % Holds all the age differences for each nSR measurement (resolution pdf)

parfor i =  1:numCores
    [nSRcounts500{i}, agediffs500{i}] = oneCoreTMRestrict(cores{i}, corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, 500);
end

%------ Run through all chosen cores with random sampling approach,
% RESTRICTION ON MINIMUM AGE DIFFERENCE = 1000

%Initiate variables
nSRcounts1000   = cell(numCores, 1); % Holds all the nSR counts (which form histogram that makes nSR pdf)
agediffs1000    = cell(numCores, 1); % Holds all the age differences for each nSR measurement (resolution pdf)

parfor i =  1:numCores
    [nSRcounts1000{i}, agediffs1000{i}] = oneCoreTMRestrict(cores{i}, corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, 1000);
end


%% Calculate transition matrices for each setup
% TM0     = TMcalculation(nSRcounts);
% TM1000  = TMcalculation(nSRcounts1000);

%% Plot histograms and nSR distributions
%Plot histograms for All Cores and 0 restriction
plotSRandResHistograms(nSRcounts, agediffs, num14cpairs, allCoresLog, 101, 'k', "All Cores")

 %Plot histograms for All Cores and 500y restriction
 plotSRandResHistograms(nSRcounts500, agediffs500, num14cpairs, allCoresLog, 101, 'k', "All Cores")

%Plot histograms for All Cores and 1000y restriction
plotSRandResHistograms(nSRcounts1000, agediffs1000, num14cpairs, allCoresLog, 101, 'k', "All Cores")

outputMetadataAndSummaryFigures(highSRhighResCoresInd, cores, lats, longs, depths, meanSR, MSI_byage, MSI_bydepth, nSRcounts, sedimentlength, num14cpairs)
