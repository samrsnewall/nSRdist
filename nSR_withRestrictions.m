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
% S - Settings to use throughout the code.

%% Add folder of necessary functions to path
addpath('Functions')

%% Create settings structure

%Set up paths to sandbox, full path to Rscript, World Atlas data on
%computer, 
S.sandboxPath  = "/Volumes/ExtDrive850X/MATLAB/nSRdist_code";
S.RscriptPath  = "/usr/local/bin/Rscript";
S.WApath       = "/Applications/PaleoDataView/WA_Foraminiferal_Isotopes_2022";
S.sheet        = "DataSheets/COPYcore40MetadataAndLin2014_2.xlsx";

%Set up file to save outputs to
stringID = "All1_RLGtrue_DS0p05_SHAK065K_Jul16";
S.dataOutputFile = "dataT_" + stringID;

%Set up core choice settings
S.minimumCoreDepth = 1000;%Minimum core depth (mbsf) in m
S.maxAtlanticLatN   = 40; %Maximum north latitude of Atlantic Cores - These not used yet
S.maxOtherLatN      = 40; %... - not used yet
S.maxAtlanticLatS   = 40; %... - not used yet
S.matOtherLatS      = 40; %... - not used yet

%Set up setting choices fo nSR calculations 
S.replicateLin2014  = false;        %See next Section for Lin2014 Set Up 
S.useLin            = true;    %Use the ages from Lin2014 database 
S.modifyLin2014Data = true;    %Whether to use dates as used by Lin2014 or to include my own modifications to remove reversals or large age gaps
S.usePF             = true;     %Use the ages I've added
S.DeltaRError       = 200;      %Error put on the Delta R (reservoir age correction)
S.c14AgeLim         = [0 50];   %Cutoffs for radiocarbon ages, in kyr
S.weighting         = "depth";  %How to weight (options are "depth", "age", or "none")
S.normWithRunMean   = true;    %Use a common meanSR to use when calculating normalised SR (false) or use the SR from each individual run (true).
S.numruns           = 1000;     %How many sets of samples to take when using the probabilistic approaches.

%Set up parameters that influence Bchron running
S.useBchron         = true;     %
S.BchronFilter      = true;     % Filter data (remove outliers and large gaps) before running Bchron or not
S.BchronFolderName  = "All1_RLGtrue_DS0p05_Jun2";    %What folder to get BchronInputs from
S.BchronOutlier     = 0.05;     %Value to input to Bchrons OutlierProbs
S.BchronReversalCriteria = 0.75; % What fraction of Bchron runs must reject an age to classify it as a reversal
S.BchronCalCurve    = "Marine20";%What calibration curve to use in Bchron
S.BchronDepthSpacing = 0.05;     %Bchron outputs the age estimates from each Bchron MCMC sample for all depths from shallowest to deepest with this interval
S.BchronReDo        = false;     %Whether to redo all Bchron regardless of whether there is already an existing Bchron run available. Options are false and true

%Set up parameters that influence Newall's Method
S.minNumberOfAges   = 4;        %Minimum number of ages a core must have (after filtering) to be used
S.reversalCriteria  = 0.75;     %What fraction of SRs between two ages must be negative to call it a reversal
S.removeLargeGaps   = true;     %Whether to manually remove large age gaps or leave them in
S.pdfMinVal         = 1e-6;     %Cutoff to reduce size of radiocarbon pdf vector to accelerate calculations
S.pdfMethod         = false;    %Whether to complete all of the pdf method or simply use it for multiply dated depths and reversals
S.useModes          = false;    %Calculate nSR distribution with the mode of each radiocarbon distribution, not sampling

%% Do I want to reproduce Lin2014 approach?
if S.replicateLin2014 == 1
    S.modifyLin2014Data = false;
    S.DeltaRError = 0;
    S.c14AgeLim = [0 55];
    S.useBchron = true;
    S.BchronCalCurve = "Marine09";
    S.BchronRedo = false;
end
%% Load Metadata of MSPF cores
%Go through sheet that has details about which dates are from planktonic
%foraminifera (to remove any benthic foraminifera dates stored in WA2022,
%or other undesirable material)
rawdata     = readtable(S.sheet);

% Select relevant information based on what subsets of data has been chosen
if S.useLin && ~S.usePF
    rawdataMSPF = rawdata(rawdata.Lin2014 == 1, :);
elseif ~S.useLin && S.usePF
    rawdataMSPF = rawdata(rawdata.WAuse == 1 & isnan(rawdata.Lin2014), :);
elseif S.useLin && S.usePF
    rawdataMSPF = rawdata(rawdata.WAuse == 1 | rawdata.Lin2014Keep == 1,:);
end

%% Create logical to choose all good cores
%Number of All Cores
numAllCores = length(rawdataMSPF.CoreName);

%Name of any problem cores (useful if wanting to exclude a single core)
problemCores = "PlacementString";%["GeoB3606-1", "GeoB7926-2", "GIK17961-2", "IOW226660-5", "KNR197-3-36GGC", "KNR197-3-53GGC"];
badLog       = ismember(string(rawdataMSPF.CoreName),problemCores);

%Restrict which cores to analyse based on metadata, such as depth or
%latitude
if  S.useLin && ~S.usePF
    latitudeRestrictionLog = ones(numAllCores, 1);
    depthRestrictionLog    = rawdataMSPF.WaterDepthM >= 1000;
else
    latitudeRestrictionLog = abs(rawdataMSPF.LatitudeDec) <= 40;
    depthRestrictionLog    = rawdataMSPF.WaterDepthM >= S.minimumCoreDepth;
end

%Therefore create logical of all good cores
restrictions       = latitudeRestrictionLog & depthRestrictionLog;
goodLog            = badLog == 0;% & restrictions == 1;

%% Create some other useful logicals
%Test a core based on it's name
 namedLog = ismember(string(rawdataMSPF.CoreName), "SHAK06-5K");   
% subsetChooser = namedLog;

%Test a subset of cores
subsetChooser = false(numAllCores,1);
subsetChooser(1:9) = true;

% %Test all cores
% subsetChooser = true(numAllCores, 1);

%% Get material data from Excel
% Get relevant metadata from Excel Spreadsheet into useful variables
chosenCoresLog = namedLog;
numCores    = sum(chosenCoresLog);
cores       = table2array(rawdataMSPF(chosenCoresLog, "CoreName")); %take list of MSPF corenames
lats        = table2array(rawdataMSPF(chosenCoresLog, "LatitudeDec"));
longs       = table2array(rawdataMSPF(chosenCoresLog, "LongitudeDec"));
depths      = table2array(rawdataMSPF(chosenCoresLog, "WaterDepthM"));
ocean       = table2array(rawdataMSPF(chosenCoresLog, "Basin"));
%distance2coast = table2array(rawdataMSPF(chosenCoresLog, "DistanceToCoast"));

[LabIDs, incDepths, excLabIDs, excDepths, dataLoc] = extract3(rawdataMSPF, chosenCoresLog, S);

%% Plot calibrated radiocarbon dates against depth
% %% Plot cores age depth models
%
%  for iPlot = (1:length(cores))
%   corePlotCal(cores{iPlot}, LabIDs{iPlot}, incDepths{iPlot}, excLabIDs{iPlot}, excDepths{iPlot}, dataLoc(iPlot), S)
% end

%% Get nSR values using Bchron (Mode and Probability)
% ------ Use Bchron outputs
bchronMode   = cell(numCores, 1);
bchronMedian = cell(numCores, 1);
bchronProb   = cell(numCores, 1);
bchronMeanSRs= NaN(numCores, 1);

%Choose if you want to get Bchron results
if S.useBchron
    %Run through Bchron method for all cores
    for i = 1:numCores
        [bchronMode{i}, bchronMedian{i}, bchronProb{i}, bchronMeanSRs(i)] = nSRBchron(cores{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, dataLoc(i), S);
    end
end
%% invSR PDF Approach
% This section runs through a quick, less-computationally-expensive
% estimator of SR to find any reversals, create scenarios, to calculate the
% meanSR and the number of dates used overall

if S.replicateLin2014
    %Code has problems when ages get too large, because they get to the
    %edge of the vector that is returned by MatCal
    S.c14AgeLim = [0 48];
end

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
end
disp("Created all scenarios")

%% Create a table that holds all useful metadata
agecoverage = sedimentlength./meanSR;
dataT = table(cores,lats, longs, depths, ocean, meanSR, ageModes,...
    sedimentlength, agecoverage, num14cpairs, MSI_byage, MSI_bydepth,...
    LabIDs, incDepths, excLabIDs, excDepths, dataLoc, core_invSRvals,...
    core_invSRprobs, corescenarios, scenario_meanSR);

 %% Plot figure showing all the ages being used
% plotAgeModes(chosenCoresLog, chosenCoresLog, ageModes, cores)
% %plotAgeModes2Subsets(ageModes, cores, highSRCoresLog, lowSRCoresLog)

%% nSR Random Sampling Approach
%This section runs a slower, more computationally-expensive estimator of SR
%to construct the SR distributions. The benefits that this method has over
%the previous method is that it can be used to create a transition matrix
%and it can be used to weight by age, and we can create distributions of
%age difference between radiocarbon pairs to test for a resolution effect.

%------ Run through all chosen cores with random sampling approach,
% RESTRICTION ON MINIMUM AGE DIFFERENCE = 0

%Initiate variables
disp("a")
nSRcounts       = cell(numCores, 1); % Holds all the nSR counts (which form histogram that makes nSR pdf)
agediffs        = cell(numCores, 1); % Holds all the age differences for each nSR measurement (resolution pdf)

for i =  1:numCores
    disp(cores{i})
    [nSRcounts{i}, agediffs{i}] = oneCoreTMRestrict(cores{i}, dataLoc(i), corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, scenario_meanSR{i}, ageModes{i}, S, 0);
end

% ------ Run through all chosen cores with random sampling approach,
% RESTRICTION ON MINIMUM AGE DIFFERENCE = 500

%Initiate variables
disp("b")
nSRcounts500   = cell(numCores, 1); % Holds all the nSR counts (which form histogram that makes nSR pdf)
agediffs500    = cell(numCores, 1); % Holds all the age differences for each nSR measurement (resolution pdf)

parfor i =  1:numCores
    [nSRcounts500{i}, agediffs500{i}] = oneCoreTMRestrict(cores{i}, dataLoc(i), corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, scenario_meanSR{i}, ageModes{i}, S, 500);
end

disp("c")

% ------ Run through all chosen cores with random sampling approach,
% RESTRICTION ON MINIMUM AGE DIFFERENCE = 1000

%Initiate variables
nSRcounts1000   = cell(numCores, 1); % Holds all the nSR counts (which form histogram that makes nSR pdf)
agediffs1000    = cell(numCores, 1); % Holds all the age differences for each nSR measurement (resolution pdf)

parfor i =  1:numCores
    [nSRcounts1000{i}, agediffs1000{i}] = oneCoreTMRestrict(cores{i}, dataLoc(i), corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, scenario_meanSR{i}, ageModes{i}, S, 1000);
end
disp("d")

% ------ Run through all chosen cores with random sampling approach,
% RESTRICTION ON MINIMUM AGE DIFFERENCE = 1500

%Initiate variables
nSRcounts1500   = cell(numCores, 1); % Holds all the nSR counts (which form histogram that makes nSR pdf)
agediffs1500    = cell(numCores, 1); % Holds all the age differences for each nSR measurement (resolution pdf)

parfor i =  1:numCores
    [nSRcounts1500{i}, agediffs1500{i}] = oneCoreTMRestrict(cores{i}, dataLoc(i), corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, scenario_meanSR{i}, ageModes{i}, S, 1500);
end
disp("e")
% ------ Run through all chosen cores with random sampling approach,
% RESTRICTION ON MINIMUM AGE DIFFERENCE = 2000

%Initiate variables
nSRcounts2000   = cell(numCores, 1); % Holds all the nSR counts (which form histogram that makes nSR pdf)
agediffs2000    = cell(numCores, 1); % Holds all the age differences for each nSR measurement (resolution pdf)

parfor i =  1:numCores
    [nSRcounts2000{i}, agediffs2000{i}] = oneCoreTMRestrict(cores{i}, dataLoc(i), corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, scenario_meanSR{i}, ageModes{i}, S, 2000);
end

%% Put all results into dataT table and save to an output folder
dataT = addvars(dataT, nSRcounts, nSRcounts500, nSRcounts1000, nSRcounts1500, nSRcounts2000, bchronMode, bchronMedian, bchronProb); 
 save(fullfile("Results", S.dataOutputFile), "dataT", "S", "rawdataMSPF")

%% Plot a few figures for a quick glance if wanted

%% ------ Define Subsets of interest from calculated metadata
% high SR and low SR cores (separated by 8cm/kyr following Lin et al., 2014)
depth1000Log    = depths > 1000;
lowSRCoresLog   = meanSR<= 8 & depth1000Log;
highSRCoresLog  = meanSR >8 & depth1000Log;
allCoresLog     = ~isnan(meanSR) & depth1000Log;

% find 10 highSR cores with most data
[~,highSRhighResCoresInd]   = maxk(num14cpairs.*highSRCoresLog, 10);
highSRhighResCoresLog       = unfind(highSRhighResCoresInd, numel(cores));

%% show metadata of a certain subset
metadataLog = allCoresLog;
outputMetadataAndSummaryFigures(allCoresLog,dataT)

%% Compare all runs vs individual runs result for a single set up
lognorm_BIGMACS = readtable("lognormal_BIGMACS.txt");                       % (use x values currently used in BIGMACS)
x = lognorm_BIGMACS.Var1';
BIGMACShist = readmatrix("BIGMACSdata/Lin2014_sedrateratio_cm_wo_NaN.txt");
BIGMACShist = BIGMACShist(:,4);
[MLN_BIGMACS, ~, gmfitBM] = fitMixLogNorm(BIGMACShist, x, 2, 5);
numruns2sample = 400;
desiredLog = highSRCoresLog;
desiredRestriction = nSRcounts500;
respectiveString = "HighSR500";
fitS.Lin2014AgeFiltering = 1;
fitS.Lin2014AgeFilter = [500 4000];
fitS.weighting = "depth"; 
fitS.chi2binN = 10;
fitS.dispChi2 = true;
fitS.mln1RunReps = 1;
fitS.mlnReps = 5;
fitS.enforceBinSizeLimits = true;

[mixLogAllruns] = plotSRandResHistograms(desiredRestriction, x, desiredLog, 3, 1, 2, 0, respectiveString,true, fitS);

figure;
hold on
plot(x, mixLogAllruns(:,2), '-r')
plot(x, MLN_BIGMACS(:,2), '--r')
plot(x, lognorm_BIGMACS.Var2, '--k')
xlim([0 5])

%% Try plotSRandResHistograms for bchronMode data
% oneCoreLog = allCoresLog;
% oneCoreLog(2:end) = 0;
[mixLogBchronMode, BchronModeHist] = plotSRandResHistograms(bchronMode, x, allCoresLog, 3, 1, 2, 0, "", 1, fitS);
[mixLogBchronMed, BchronMedHist] = plotSRandResHistograms(bchronMedian, x, allCoresLog, 3, 1, 2, 0, "", 1, fitS);
[mixLogBchronProb, BchronProbHist] = plotSRandResHistograms(bchronProb, x, allCoresLog, 3, 1, 2, 0, "", 1, fitS);

figure;
subplot(2,1,1)
hold on
histogram(BIGMACShist, 'Normalization','pdf', 'FaceColor', [0.5, 0.5, 0.5])
plot(x, lognorm_BIGMACS.Var2, '-k', 'LineWidth',2)
xlim([0 6])
xlabel("nSR")
title("Original BIGMACS")
subplot(2,1,2)
hold on
histogram(BchronModeHist, 'Normalization','pdf')
plot(x, lognorm_BIGMACS.Var2, '-k', 'LineWidth',2)
plot(x, mixLogBchronMode(:,2), '-r', 'LineWidth',2)
xlim([0 6])
xlabel("nSR")
title("Repeated by Newall")

figure;
hold on
histogram(BchronProbHist, 'Normalization', 'pdf')
plot(x, mixLogBchronProb(:,2))
xlim([0 10])
xlabel("nSR")

figure;
hold on
plot(x, mixLogBchronMode(:,2), '-r', 'LineWidth',2)
plot(x, mixLogBchronProb(:,2), '-b', 'LineWidth',2)
plot(x, lognorm_BIGMACS.Var2, '-k', 'LineWidth',2)
xlim([0 6])