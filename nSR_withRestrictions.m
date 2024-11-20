%% Explanation
%This script reads in the metadata xlsx sheet, accesses cores through the
%World Atlas, and then runs through the scenario creation and nSR counting
%functions. The metadata and scenario creation has the effect of choosing
%useful data, by removing reversals or creating multiple potential
%radiocarbon profiles to avoid reversals, removing data that is spaced too far apart, and deals with doubly dated
%depths. The nSR counting function takes each core and randomly samples
%many possible age samples from each dated depth to get estimates of nSR. It
%also allows a restriction to be imposed on the minimum age difference estimated
%between different radiocarbon dates.

%%%%INPUTS
% sheet - which metadata xlsx sheet to use
%
% problemCores - cores to avoid using because of poor data or excessive
% computational expense
% 

%% Add folder of necessary functions to path
addpath('Functions')

%% Load Metadata of MSPF cores
%Check which cores have MSPF (monospecific planktonic foram) dates
%sheet = "DataSheets/COPYcorechoices_MSPF_highRes.xlsx";
%sheet = "DataSheets/COPYcore40Metadata.xlsx";
sheet = "DataSheets/COPYcore40MetadataAndLin2014.xlsx";
LinOnly =0;
PFOnly = 0;
rawdata     = readtable(sheet);
if ~contains(sheet, "COPYcore40Metadata")
    rawdataMSPF = rawdata(rawdata.MSPF == 1,:);
elseif ~contains(sheet, "AndLin2014")
    rawdataMSPF = rawdata(ismember(rawdata.UseChoice, "PF"), :);
else
    if LinOnly == 1
        rawdataMSPF = rawdata(rawdata.Lin2014 == 1, :);
    elseif PFOnly == 1
        rawdataMSPF = rawdata(rawdata.WAuse == 1);
    else
        rawdataMSPF = rawdata(rawdata.WAuse == 1 | rawdata.Lin2014Keep == 1,:);
    end
end

%% Create settings structure
S.minNumberOfAges   = 4;        %Minimum number of ages a core must have (after filtering) to be used
S.DeltaRError       = 0;      %Error put on the Delta R (reservoir age correction)
S.reversalCriteria  = 0.75;     %What fraction of SRs between two ages must be negative to call it a reversal
S.removeLargeGaps   = true;     %Whether to manually remove large age gaps or leave them in
S.modifyLin2014Data = true;     %Whether to use dates as used by Lin2014 or to include my own modifications to remove reversals or large age gaps
S.c14AgeLim         = [0 42];   %Cutoffs for radiocarbon ages, in kyr %Choose [1 42] normally
S.pdfMinVal         = 1e-6;     %Cutoff to reduce size of radiocarbon pdf vector to accelerate calculations
S.pdfMethod         = false;    %Whether to complete all of the pdf method or simply use it for multiply dated depths and reversals
S.weighting         = true;     %Whether or not to weight by depth (= false means no weighting)
S.normWithRunMean   = false;    %Use a common meanSR to use when calculating normalised SR (false) or use the SR from each individual run (true).
S.useModes          = false;     %Calculate nSR distribution with the mode of each radiocarbon distribution, not sampling

%% Choose cores to look at
%------ Index Good Cores
%Create index for all cores
numAllCores = length(rawdataMSPF.CoreName);
allcores    = 1:numAllCores;

%Create logical of cores without those that have been manually removed
%(goodLog)
reversalDenseCores = [ "H214", "SO75_3_26KL", "KNR159-5-36GGC", "MD02-2550", "GIK17940-2", "MD07-3076"];
if LinOnly == 1
    reversalDenseCores = "PLACEMENTSTRING";
end
problemCores       = []; %reversals are mixed with duplicated depths, haven't figured out how to handle this yet
badLog             = ismember(string(rawdataMSPF.CoreName),[reversalDenseCores, problemCores]);
badIndexes         = allcores(badLog);
goodLog            = badLog == 0;
goodIndexes = allcores(~badLog);

%Create some other useful logical, to test all cores, only 1 core (defined by name), or a
%subset of cores (defined by index)
namedLog           = ismember(string(rawdataMSPF.CoreName), "DSDP594");
allLog             = true(length(goodLog), 1);
subsetChooser      = logical(zeros(numAllCores, 1)); %#ok<LOGL>
subsetChooser(1:end)  = 1; 
subsetChooser = subsetChooser == 1 & goodLog == 1;

%% Get material data from Excel
% Get info from Excel Spreadsheet
chosenCoresLog = goodLog;
numCores    = sum(chosenCoresLog);
cores       = table2array(rawdataMSPF(chosenCoresLog, "CoreName")); %take list of MSPF corenames
lats        = table2array(rawdataMSPF(chosenCoresLog, "LatitudeDec"));
longs       = table2array(rawdataMSPF(chosenCoresLog, "LongitudeDec"));
depths      = table2array(rawdataMSPF(chosenCoresLog, "WaterDepthM"));
ocean       = table2array(rawdataMSPF(chosenCoresLog, "Basin"));
%distance2coast = table2array(rawdataMSPF(chosenCoresLog, "DistanceToCoast"));
if ~contains(sheet, "COPYcore40Metadata")
[LabIDs, incDepths, excLabIDs, excDepths] = extract1(rawdataMSPF, chosenCoresLog);
dataLoc = strings(numCores,1);
dataLoc(:) = "WA";

elseif ~contains(sheet, "AndLin2014")
    [LabIDs, incDepths, excLabIDs, excDepths] = extract2(rawdataMSPF, chosenCoresLog, S);
    dataLoc = strings(numCores,1);
dataLoc(:) = "WA";
else
    [LabIDs, incDepths, excLabIDs, excDepths, dataLoc] = extract3(rawdataMSPF, chosenCoresLog, LinOnly, S);
end


% %% Plot cores age depth models
% %------ Plot radiocarbon data from a single core (for exploratory use)
% % for iPlot = 1
%  for iPlot = (1:length(cores))
% 
%   corePlotCal(cores{iPlot}, LabIDs{iPlot}, incDepths{iPlot}, excLabIDs{iPlot}, excDepths{iPlot}, S)
% end

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
ageModes        = cell(numCores, 1);
transprobs_cores= nan(3,3,numCores);
corescenarios   = cell(numCores,1);
newlabels       = cell(numCores,1);
numreversals    = nan(numCores,1);
scenario_meanSR = cell(numCores, 1);

%% ----- Apply invSR with loop
%Calculate SR distribution for each core, as well as meanSR and other
%useful information
for i = 1:numCores
    disp("Working on " + cores{i})
    [core_invSRvals{i}, core_invSRprobs{i}, meanSR(i), MSI_byage(i),...
        MSI_bydepth(i), sedimentlength(i), num14cpairs(i), ageModes{i},...
        corescenarios{i}, newlabels{i}, numreversals(i), scenario_meanSR{i}]...
        = oneCoreSRpdf(cores{i}, dataLoc(i), LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, S, 0);
    disp("Done with " + cores{i})
end
disp("Created all scenarios")

%% Create a table that holds all useful metadata
agecoverage = sedimentlength./meanSR;
dataT = table(cores,lats, longs, depths, ocean, meanSR, ageModes, sedimentlength, agecoverage, num14cpairs, MSI_byage, MSI_bydepth, LabIDs, incDepths, excLabIDs, excDepths, core_invSRvals, core_invSRprobs, corescenarios, scenario_meanSR);

%% ------ Define Subsets of interest
% high SR and low SR cores (separated by 8cm/kyr following Lin et al., 2014)
depth1000Log = depths > 1000;
lowSRCoresLog   = meanSR<= 8 & depth1000Log;
highSRCoresLog  = meanSR >8 & depth1000Log;
allCoresLog     = ~isnan(meanSR) & depth1000Log;

% find 10 highSR cores with most data
[~,highSRhighResCoresInd]   = maxk(num14cpairs.*highSRCoresLog, 10);
highSRhighResCoresLog       = unfind(highSRhighResCoresInd, numel(cores));

%% plot figure showing all the ages being used
%plotAgeModes(depth1000Log, ageModes, cores)
%plotAgeModes2Subsets(ageModes, cores, highSRCoresLog, lowSRCoresLog)

%% nSR Random Sampling Approach
%This section runs a slower, more computationally-expensive estimator of SR
%to construct the SR distributions. The benefits that this method has over
%the previous method is that it can be used to create a transition matrix
%and it can be used to weight by age, and we can create distributions of
%age difference between radiocarbon pairs to test for a resolution effect.

%------ Run through all chosen cores with random sampling approach,
% RESTRICTION ON MINIMUM AGE DIFFERENCE = 0

% for i = 1:numCores
%     corescenarios{i} = labels{i}
% end

%Initiate variables
disp("a")
nSRcounts       = cell(numCores, 1); % Holds all the nSR counts (which form histogram that makes nSR pdf)
agediffs        = cell(numCores, 1); % Holds all the age differences for each nSR measurement (resolution pdf)

for i =  1:numCores
    [nSRcounts{i}, agediffs{i}] = oneCoreTMRestrict(cores{i}, dataLoc(i), corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, scenario_meanSR{i}, ageModes{i}, S, 0);
end

% ------ Run through all chosen cores with random sampling approach,
% RESTRICTION ON MINIMUM AGE DIFFERENCE = 500

%Initiate variables
disp("b")
nSRcounts500   = cell(numCores, 1); % Holds all the nSR counts (which form histogram that makes nSR pdf)
agediffs500    = cell(numCores, 1); % Holds all the age differences for each nSR measurement (resolution pdf)

for i =  1:numCores
    [nSRcounts500{i}, agediffs500{i}] = oneCoreTMRestrict(cores{i}, dataLoc(i), corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, scenario_meanSR{i}, S, 500);
end

disp("c")

% ------ Run through all chosen cores with random sampling approach,
% RESTRICTION ON MINIMUM AGE DIFFERENCE = 1000

%Initiate variables
nSRcounts1000   = cell(numCores, 1); % Holds all the nSR counts (which form histogram that makes nSR pdf)
agediffs1000    = cell(numCores, 1); % Holds all the age differences for each nSR measurement (resolution pdf)

parfor i =  1:numCores
    [nSRcounts1000{i}, agediffs1000{i}] = oneCoreTMRestrict(cores{i}, dataLoc(i), corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, scenario_meanSR{i}, S, 1000);
end
disp("d")

% ------ Run through all chosen cores with random sampling approach,
% RESTRICTION ON MINIMUM AGE DIFFERENCE = 1000

%Initiate variables
nSRcounts1500   = cell(numCores, 1); % Holds all the nSR counts (which form histogram that makes nSR pdf)
agediffs1500    = cell(numCores, 1); % Holds all the age differences for each nSR measurement (resolution pdf)

parfor i =  1:numCores
    [nSRcounts1500{i}, agediffs1500{i}] = oneCoreTMRestrict(cores{i}, dataLoc(i), corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, scenario_meanSR{i}, S, 1500);
end
disp("e")
% ------ Run through all chosen cores with random sampling approach,
% RESTRICTION ON MINIMUM AGE DIFFERENCE = 1000

%Initiate variables
nSRcounts2000   = cell(numCores, 1); % Holds all the nSR counts (which form histogram that makes nSR pdf)
agediffs2000    = cell(numCores, 1); % Holds all the age differences for each nSR measurement (resolution pdf)

for i =  1:numCores
    [nSRcounts2000{i}, agediffs2000{i}] = oneCoreTMRestrict(cores{i}, dataLoc(i), corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, scenario_meanSR{i}, S, 2000);
end

dataT = addvars(dataT, nSRcounts, nSRcounts500, nSRcounts1000, nSRcounts1500, nSRcounts2000); 

metadataLog = allCoresLog;
outputMetadataAndSummaryFigures(allCoresLog,dataT)

%save("Results/dataT_PFAndLin2014_Nov19S1.mat", "dataT", "S")

%% Calculate transition matrices for each setup
% [~,~,TM1000]  = TMcalculation(nSRcounts1000, highSRCoresLog);

%x = 0.01:0.01:15; plotSRandResHistograms(dataT.nSRcounts500, x, ones(size(dataT.nSRcounts500)), true, 3,1,2,0,"",true)

%% Plot a few figures for a quick glance if wanted

%% Compare all runs vs individual runs result for a single set up
lognorm_BIGMACS = readtable("lognormal_BIGMACS.txt");                       % (use x values currently used in BIGMACS)
x = lognorm_BIGMACS.Var1';
numruns2sample = 400;
desiredLog = highSRCoresLog;
desiredRestriction = nSRcounts;
respectiveString = "HighSR500";
fitS.Lin2014AgeFiltering = 1;

[mixLogAllruns] = plotSRandResHistograms(desiredRestriction, x, desiredLog, true, 3, 1, 2, 0, respectiveString,true, fitS);

figure;
hold on
plot(x, mixLogAllruns(:,2), '-r')
plot(x, lognorm_BIGMACS.Var2)
xlim([0 5])
