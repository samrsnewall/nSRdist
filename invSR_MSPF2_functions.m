%% ReadMe
%This is the main script for creating a distribution of normalised
%sedimentation rate from radiocarbon dated ocean sediment cores. The
%radiocarbon data is taken from the World Atlas of late Quaternary
%Foraminifera Oxygen and Carbon Isotope Ratios (Mulitza et al., 2022) -
%henceforth called WA2022. The data is filtered using the excel file
%COPYcorechoices_MSPF.xlsx, which contains the information needed to choose
%only the monospecific planktonic foraminifera - henceforth MSPF dates from
%each core. This excel file also allows for some dates to be removed "by
%hand" to speed up computational process.

%A few of the important choices made are:
% - If >3/4 of possible pathways between pair of radiocarbon dates is negative, it is classified as a reversal
% - The \DeltaR is 0±200yr (change this in multiMatcal.m function)
% - Cores are split by SR at 8cm/kyr
% - An extra subset is made by picking 10 highSR cores with most data

% Add Functions Folder to Path
addpath('Functions')

%% Load Metadata of MSPF cores
%Check which cores have MSPF (monospecific planktonic foram) dates
data = readtable("COPYcorechoices_MSPF_highRes2.xlsx"); %read all metadata
dataMSPF = data(data.MSPF == 1,:);

%------ Index Good Cores
%Create index for all cores
numAllCores = length(dataMSPF.CoreName);
allcores = 1:numAllCores;
%Create index of cores that have many reversals (determined by manual
%inspection)
reversalDenseCores = ["GeoB1711-4", "H214", "SO75_3_26KL", "KNR159-5-36GGC"];
%problemCores = ["SU81-18", "MD88-770"]; % These all have too many reversals 
problemCores = [];
badLog = contains(string(dataMSPF.CoreName),[reversalDenseCores, problemCores]);
namedLog = contains(string(dataMSPF.CoreName), "GEOFARKF13");
goodLog = badLog == 0;
allLog = true(length(goodLog), 1);
badIndexes = allcores(badLog);

%Hence create index to use only cores that passed manual inspection
goodIndexes = allcores(~badLog);

%------- Create a subset of cores to look at if interested
%Decide which subset of cores to look at
subsetChooser = logical(zeros(numAllCores, 1)); %#ok<LOGL>
subsetChooser(1:5) = 1; 
subsetChooser(badLog) = 0;

%------- Take desired data into arrays
chosenCoresLog = goodLog;
cores = table2array(dataMSPF(chosenCoresLog, "CoreName")); %take list of MSPF corenames
lats = table2array(dataMSPF(chosenCoresLog, "LatitudeDec"));
longs = table2array(dataMSPF(chosenCoresLog,"LongitudeDec"));
depths = table2array(dataMSPF(chosenCoresLog, "WaterDepthM"));
LabIDs = table2cell(dataMSPF(chosenCoresLog, "LabIDs")); %take list of LabIDs relating to MSPF dates of each core
incDepths = table2cell(dataMSPF(chosenCoresLog, "IncludeDepths")); % take list of depths (useful if no labels)
excLabIDs = table2cell(dataMSPF(chosenCoresLog, "excludeLabIDs")); %take list of manually removed dates for each core
excDepths = table2cell(dataMSPF(chosenCoresLog, "excludeDepth")); %take list of manually removed dates for each core (useful if no labels)
numCores = sum(chosenCoresLog);

%------ Plot radiocarbon data from a single core (for exploratory use)
%for iPlot = find(chosenCoresLog)'
for iPlot = (1:length(incDepths))
  corePlot(cores{iPlot}, LabIDs{iPlot}, incDepths{iPlot}, excLabIDs{iPlot}, excDepths{iPlot})
end

%% invSR PDF Approach
% This section runs through a quick, less-computationally-expensive
% estimator of SR to find any reversals, create scenarios, to calculate the
% meanSR and the number of dates used overall

%------ Initialise variables to hold this information
core_invSRvals = cell(numCores,1);
core_invSRprobs = cell(numCores,1);
meanSR = nan(numCores,1);
MSI_byage = nan(numCores,1);
MSI_bydepth = nan(numCores,1);
sedimentlength = nan(numCores,1);
num14cpairs = nan(numCores,1);
transprobs_cores= nan(3,3,numCores);
corescenarios = cell(numCores,1);
newlabels = cell(numCores,1);
numreversals = nan(numCores,1);

%----- Apply invSR with loop
%Calculate SR distribution for each core, as well as meanSR and other
%useful information

for i = 1:numCores
    disp(cores{i})
    [core_invSRvals{i}, core_invSRprobs{i}, meanSR(i), MSI_byage(i), MSI_bydepth(i), sedimentlength(i), num14cpairs(i), corescenarios{i}, newlabels{i}, numreversals(i)] = oneCoreSRpdf(cores{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, 0);
end

%------ Define Subsets of interest
% high SR and low SR cores (separated by 8cm/kyr following Lin et al., 2014)
lowSRCoresLog = meanSR<= 8;
highSRCoresLog = meanSR >8;
allCoresLog = ~isnan(meanSR);

% find 10 highest resolution
% MSI_byageTOOL = MSI_byage.*highSRCoresLog; MSI_byageTOOL(isnan(MSI_byage) | MSI_byageTOOL == 0) = 9e9;
% [~,highSRhighResCoresInd] = maxk(-MSI_byageTOOL, 10);

% find 10 highSR cores with most data
[~,highSRhighResCoresInd] = maxk(num14cpairs.*highSRCoresLog, 10);
highSRhighResCoresLog = unfind(highSRhighResCoresInd, numel(cores));

%------ Combine the invSR pdfs from each core
[master_invSRvals, allcore_invSRprobs] = combinepdfs(core_invSRvals, core_invSRprobs, sedimentlength);

%------ Convert to Sed Rate
%Change inverseSR distribution to SR distribution
[SRvals_interp, SRvalsprob_norm] = invSRtoSR(master_invSRvals, allcore_invSRprobs);

%------ Plot the SR distribution from PDF method and compare with BIGMACS values

%Load normalised SR data from BIGMACS
normSR_BIGMACS = load("../BIGMACSdata/Lin2014_sedrateratio_cm_wo_NaN.txt");
%Load lognormal file
lognormdata_BIGMACS = load("../BIGMACSdata/lognormal.txt");

figure;
subplot(2,1,1)
yyaxis left
ylabel("Depth Weighted Counts")
histogram(normSR_BIGMACS, "FaceColor", [0.8 0.8 0.8])
hold on
yyaxis right
plot(lognormdata_BIGMACS(:,1), lognormdata_BIGMACS(:,2))
ylim([0 1.1])
ylabel("Probability")
xlabel("Sed Rate Ratio")
xlim([0 6])
legend("BIGMACS Histogram")

subplot(2,1,2)
plot(SRvals_interp, SRvalsprob_norm)
xlim([0 6])
ylabel("Probability")
xlabel("Sed Rate Ratio")
legend("MSPF2 PDF")


%% invSR Random Sampling Approach
%This section runs a slower, more computationally-expensive estimator of SR
%to construct the SR distributions. The benefits that this method has over
%the previous method is that it can be used to create a transition matrix
%and it can be used to weight by age, and we can create distributions of
%age difference between radiocarbon pairs to test for a resolution effect.

%----- Initiate variables
nSRcounts = cell(numCores,1); % Holds all the nSR counts (which form histogram that makes nSR pdf)
agediffs = cell(numCores, 1); % Holds all the age differences for each nSR measurement (resolution pdf)
coreTM = nan(3,3,numCores);
cores_transnums = nan(3,3,numCores); %Counts of each type of transition
cores_CSE2x = nan(3,1,numCores); %Counts of the number of transitions from each starting state

%------ Run through all chosen cores with random sampling approach
calcTM = false; %do you want TM code to be run? (False if no, True if yes)
for i =  1:numCores
    disp(cores{i})
    [cores_transnums(:,:,i), cores_CSE2x(:,:,i), nSRcounts{i}, agediffs{i}] = oneCoreTM(cores{i}, corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, calcTM);
end

if calcTM ~= false
    %----- Calculate all core TM
    %Calculate the number of transitions across all the cores
    allcores_transnums = sum(cores_transnums, 3, 'omitmissing');
    allcores_CSE2x = sum(cores_CSE2x, 3, 'omitmissing');
    allcores_CSE2x_ratios = allcores_CSE2x'./(sum(allcores_CSE2x));

    %----- Construct the transition matrix
    allTM = [allcores_CSE2x_ratios;allcores_transnums./allcores_CSE2x];
end

%% Plot histograms of random sample counts
%Plot histograms for All Cores
plotSRandResHistograms(nSRcounts, agediffs, num14cpairs, allCoresLog, 101, 'k', "All Cores")

%Plot histograms for High SR subset
if sum(highSRCoresLog > 0)
plotSRandResHistograms(nSRcounts, agediffs, num14cpairs, highSRCoresLog, 101, 'r', "High SR")
end

%Plot histograms for Low SR Subset
if sum(lowSRCoresLog >0)
plotSRandResHistograms(nSRcounts, agediffs, num14cpairs, lowSRCoresLog, 102, 'b', "Low SR")
end

%Plot 
if sum(highSRCoresLog) >10
    plotSRandResHistograms(nSRcounts, agediffs, num14cpairs, highSRhighResCoresLog, 103, 'k', "10 Highest Res")
end

%% Plot output metadata and figures

%----- Find index of all cores used in the analysis
coresUsedLog = ~cellfun('isempty', core_invSRvals);
coresUsedInd = find(coresUsedLog);

%----- Get output metadata and summary figures for all cores used
outputMetadataAndSummaryFigures(coresUsedInd, cores, lats, longs, depths, meanSR, MSI_byage, MSI_bydepth, nSRcounts, sedimentlength, num14cpairs)



