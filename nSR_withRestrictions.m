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
% problemCores - cores to avoid using because of poor data or excessive
% computational expense
% S - Settings to use throughout the code
% 

%% Add folder of necessary functions to path
addpath('Functions')

%% Create settings structure

%Set up paths to sandbox, full path to Rscript, and World Atlas data on
%computer
S.sandboxPath       = "/Volumes/ExtDrive850X/MATLAB/nSRdist_code";
S.RscriptPath       = "/usr/local/bin/Rscript";
S.WApath            = "/Applications/PaleoDataView/WA_Foraminiferal_Isotopes_2022";

%Set up setting choices for nSR calculations
S.minNumberOfAges   = 4;        %Minimum number of ages a core must have (after filtering) to be used
S.DeltaRError       = 200;      %Error put on the Delta R (reservoir age correction)
S.reversalCriteria  = 0.75;     %What fraction of SRs between two ages must be negative to call it a reversal
S.removeLargeGaps   = false;    %Whether to manually remove large age gaps or leave them in
S.modifyLin2014Data = false;    %Whether to use dates as used by Lin2014 or to include my own modifications to remove reversals or large age gaps
S.c14AgeLim         = [0 45];   %Cutoffs for radiocarbon ages, in kyr %Choose [1 42] normally
S.pdfMinVal         = 1e-6;     %Cutoff to reduce size of radiocarbon pdf vector to accelerate calculations
S.pdfMethod         = false;    %Whether to complete all of the pdf method or simply use it for multiply dated depths and reversals
S.weighting         = "depth";  %How to weight (options are "depth", "age", or "none")
S.normWithRunMean   = false;    %Use a common meanSR to use when calculating normalised SR (false) or use the SR from each individual run (true).
S.useModes          = false;    %Calculate nSR distribution with the mode of each radiocarbon distribution, not sampling
S.useBchron         = true;     %
S.useLin            = true;     %Use the ages from Lin2014 database 
S.usePF             = true;    %Use the ages I've added
S.BchronFolderName  = '-';      %What folder to get BchronInputs from
S.BchronOutlier     = 0.05;     %Value to input to Bchrons OutlierProbs
S.BchronCalCurve    ="Marine09";%What calibration curve to use in Bchron
S.BchronReDo        = true;     %Whether to redo all Bchron regardless of whether there is already an existing Bchron run available. Options are false and true

%% Do I want to replicate Lin2014 approach?

S.replicateLin2014 = 1;

if S.replicateLin2014 == 1
    S.modifyLin2014Data = false;
    S.DeltaRError = 0;
    S.c14AgeLim = [0 48]; %Code has problems when age gets too large, gets to edge of calibration curve
    S.useBchron = true;
    S.BchronCalCurve = "Marine09";
end


%% Load Metadata of MSPF cores
%Check which cores have MSPF (monospecific planktonic foram) dates
sheet = "DataSheets/COPYcore40MetadataAndLin2014.xlsx";
rawdata     = readtable(sheet);

% Select relevant information based on what subsets of data has been chosen
if S.useLin && ~S.usePF
    rawdataMSPF = rawdata(rawdata.Lin2014 == 1, :);
elseif ~S.useLin && S.usePF
    rawdataMSPF = rawdata(rawdata.WAuse == 1, :);
elseif S.useLin && S.usePF
    rawdataMSPF = rawdata(rawdata.WAuse == 1 | rawdata.Lin2014Keep == 1,:);
end

%% Create logical to choose all good cores
%Number of All Cores
numAllCores = length(rawdataMSPF.CoreName);

%Name of any problem cores (useful if wanting to exclude a single core)
problemCores = "PLACEMENTSTRING";
badLog             = ismember(string(rawdataMSPF.CoreName),problemCores);

%Restrict which cores to analyse based on metadata, such as depth or
%latitude
if  S.useLin && ~S.usePF
    latitudeRestrictionLog = ones(numAllCores, 1);
    depthRestrictionLog    = rawdataMSPF.WaterDepthM >= 1000;
else
    latitudeRestrictionLog = abs(rawdataMSPF.LatitudeDec) <= 40;
    depthRestrictionLog    = rawdataMSPF.WaterDepthM >= 1000;
end

%Therefore create logical of all good cores
restrictions       = latitudeRestrictionLog & depthRestrictionLog;
goodLog            = badLog == 0 & restrictions == 1;


%% Create some other useful logicals
%Test a core based on it's name
namedLog = ismember(string(rawdataMSPF.CoreName), "H214");   

%Test a subset of cores
subsetChooser = false(numAllCores, 1);
subsetChooser(1:end) = 1; 
subsetChooser = subsetChooser == 1 & goodLog == 1;

%% Get material data from Excel
% Get relevant metadata from Excel Spreadsheet into useful variables
chosenCoresLog = goodLog;
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
%   corePlotCal(cores{iPlot}, LabIDs{iPlot}, incDepths{iPlot}, excLabIDs{iPlot}, excDepths{iPlot}, S)
% end

%% Get nSR values using Bchron (Mode and Probability)
% ------ Use Bchron outputs
bchronMode   = cell(numCores, 1);
bchronMedian = cell(numCores, 1);
bchronProb   = cell(numCores, 1);

for i = 1:numCores
[bchronMode{i}, bchronMedian{i}, bchronProb{i}] = nSRBchron(cores{i}, dataLoc(i), S);
end

%% Get transition matrix using Bchron Mode
 [~,~,TMBchronMode]  = TMcalculation(bchronMode, true(numCores,1), S);

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
dataT = table(cores,lats, longs, depths, ocean, meanSR, ageModes, sedimentlength, agecoverage, num14cpairs, MSI_byage, MSI_bydepth, LabIDs, incDepths, excLabIDs, excDepths, core_invSRvals, core_invSRprobs, corescenarios, scenario_meanSR);

%% ------ Define Subsets of interest from calculated metadata
% high SR and low SR cores (separated by 8cm/kyr following Lin et al., 2014)
depth1000Log = depths > 1000;
lowSRCoresLog   = meanSR<= 8 & depth1000Log;
highSRCoresLog  = meanSR >8 & depth1000Log;
allCoresLog     = ~isnan(meanSR) & depth1000Log;

% find 10 highSR cores with most data
[~,highSRhighResCoresInd]   = maxk(num14cpairs.*highSRCoresLog, 10);
highSRhighResCoresLog       = unfind(highSRhighResCoresInd, numel(cores));

 %% Plot figure showing all the ages being used
% %plotAgeModes(depth1000Log, ageModes, cores)
%plotAgeModes2Subsets(ageModes, cores, highSRCoresLog, lowSRCoresLog)

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

%%

%Initiate variables
disp("b")
nSRcounts500   = cell(numCores, 1); % Holds all the nSR counts (which form histogram that makes nSR pdf)
agediffs500    = cell(numCores, 1); % Holds all the age differences for each nSR measurement (resolution pdf)

for i =  1:numCores
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

%%

dataT = addvars(dataT, nSRcounts, nSRcounts500, nSRcounts1000, nSRcounts1500, nSRcounts2000, bchronMode, bchronMedian, bchronProb); 

metadataLog = allCoresLog;
outputMetadataAndSummaryFigures(allCoresLog,dataT)

save("Results/dataT_LinandPF_LinMethod_Dec10.mat", "dataT", "S")

%% Calculate transition matrices for each setup


%x = 0.01:0.01:15; plotSRandResHistograms(dataT.nSRcounts500, x, ones(size(dataT.nSRcounts500)), true, 3,1,2,0,"",true)

%% Plot a few figures for a quick glance if wanted

%% Compare all runs vs individual runs result for a single set up
lognorm_BIGMACS = readtable("lognormal_BIGMACS.txt");                       % (use x values currently used in BIGMACS)
x = lognorm_BIGMACS.Var1';
BIGMACShist = readmatrix("BIGMACSdata/Lin2014_sedrateratio_cm_wo_NaN.txt");
BIGMACShist = BIGMACShist(:,4);
[MLN_BIGMACS, ~, gmfitBM] = fitMixLogNorm(BIGMACShist, x, 2, 0, 5);
numruns2sample = 400;
desiredLog = highSRCoresLog;
desiredRestriction = nSRcounts500;
respectiveString = "HighSR500";
fitS.Lin2014AgeFiltering = 1;
fitS.weighting = "depth";
fitS.chi2binN = 10;

[mixLogAllruns] = plotSRandResHistograms(desiredRestriction, x, desiredLog, true, 3, 1, 2, 0, respectiveString,true, fitS);

figure;
hold on
%plot(x, mixLogAllruns(:,2), '-r')
plot(x, MLN_BIGMACS(:,2), '--r')
plot(x, lognorm_BIGMACS.Var2, '--k')
xlim([0 5])

%% Try plotSRandResHistograms for bchronMode data
% oneCoreLog = allCoresLog;
% oneCoreLog(2:end) = 0;
[mixLogBchronMode, BchronModeHist] = plotSRandResHistograms(bchronMode, x, allCoresLog, 1, 3, 1, 2, 0, "", 1, fitS);
[mixLogBchronMed, BchronMedHist] = plotSRandResHistograms(bchronMedian, x, allCoresLog, 1, 3, 1, 2, 0, "", 1, fitS);
[mixLogBchronProb, BchronProbHist] = plotSRandResHistograms(bchronProb, x, allCoresLog, 1, 3, 1, 2, 0, "", 1, fitS);

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