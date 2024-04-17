%% ReadMe
%We use corechoices_MSPF.xlsx, which has metadata for each core, including
%the labIDs of the samples that are MSPF of the same species, and including
%information on which data we believe should be removed (from visual
%inspection).

%% Add Functions Folder to Path
addpath('Functions')

%% Load Metadata of MSPF cores
%Check which cores have MSPF (monospecific planktonic foram) dates
data = readtable("../../Data/corechoices_MSPF.xlsx"); %read all metadata
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
subsetchooser(50:80) = 1;
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

%% invSR PDF Approach
% This section runs through a quicker, less-computationally-expensive
% estimator of SR to find any reversals (and thence create scenarios to
% avoid them), to calculate the meanSR, return number of dates used. 

% Initialise variables to hold this information
core_invSRvals = cell(1,numCores);
core_invSRprobs = cell(1,numCores);
meanSR = nan(1,numCores);
res_byage = nan(1,numCores);
res_bydepth = nan(1,numCores);
sedimentlength = nan(1,numCores);
num14cpairs = nan(1,numCores);
transprobs_cores= nan(3,3,numCores);
corescenarios = cell(1,numCores);
newlabels = cell(1,numCores);
numreversals = nan(1,numCores);

%Calculate SR distribution for each core, as well as meanSR and other
%useful information
parfor i = 1:numCores
    disp(cores{i})
    [core_invSRvals{i}, core_invSRprobs{i}, meanSR(i), res_byage(i), res_bydepth(i), sedimentlength(i), num14cpairs(i), corescenarios{i}, newlabels{i}, numreversals(i)] = oneCoreSRpdf(cores{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, 0);
end

%% Find high SR and low SR cores (separated by 8cm/kyr following Lin et al., 2014)
lowSRCoresInd = meanSR<= 8;
highSRCoresInd = meanSR >8;

% find top 10 highest resolution cores
[~,highRes10CoresInd] = maxk(-res_byage, 10);
cores(highRes10CoresInd)

%% Combine the invSR pdfs from each core
[master_invSRvals, allcore_invSRprobs] = combinepdfs(core_invSRvals, core_invSRprobs, sedimentlength);

%% Convert to Sed Rate
%Change inverseSR distribution to SR distribution
[SRvals_interp, SRvalsprob_norm] = invSRtoSR(master_invSRvals, allcore_invSRprobs);

%% Plot the SR distribution from PDF method and compare with BIGMACS values
%Load normalised SR data from BIGMACS
normSR_BIGMACS = load("../../Data/Lin2014_sedrateratio_cm_wo_NaN.txt");
%Load lognormal file
lognormdata_BIGMACS = load("../lognormal.txt");

figure(22)
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

%Initiate variables
nSRcounts = cell(numCores,1); % Holds all the nSR counts (which form histogram that makes nSR pdf)
agediffs = cell(numCores, 1); % Holds all the age differences for each nSR measurement (resolution pdf)
coreTM = nan(3,3,numCores);
cores_transnums = nan(3,3,numCores); %Counts of each type of transition
cores_CSE2x = nan(3,1,numCores); %Counts of the number of transitions from each starting state

%Run through all chosen cores with random sampling approach
parfor i =  1:numCores
    disp(cores{i})
    [cores_transnums(:,:,i), cores_CSE2x(:,:,i), nSRcounts{i}, agediffs{i}] = oneCoreTM(cores{i}, corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i});
end

%% Calculate all core TM
%Calculate the number of transitions across all the cores
allcores_transnums = sum(cores_transnums, 3, 'omitmissing');
allcores_CSE2x = sum(cores_CSE2x, 3, 'omitmissing');
allcores_CSE2x_ratios = allcores_CSE2x'./(sum(allcores_CSE2x));

%Construct the transition matrix
allTM = [allcores_CSE2x_ratios;allcores_transnums./allcores_CSE2x];

%% Plot the histogram of random sample counts
plotSRandResHistograms(nSRcounts, agediffs, num14cpairs, highSRCoresInd, 101, 'k', "High SR")
plotSRandResHistograms(nSRcounts, agediffs, num14cpairs, lowSRCoresInd, 102, 'k', "Low SR")
plotSRandResHistograms(nSRcounts, agediffs, num14cpairs, highRes10CoresInd, 103, 'k', "High Res 10")
figure(102)
subplot(1,2,2)
ylim([0 100])

%% Plot Maps of data in MATLAB
%Find index of all cores used in the analysis
ind2 = NaN(1,numCores);
for i = 1:numCores
    if isempty(core_invSRvals{i})
        ind2(i) = 0;
    else
        ind2(i) =i;
    end
end
ind3 = ind2(ind2 ~= 0);

%Find corenames, lats, longs, depths of cores included
core_inc = cores(ind3);
lat_inc = lats(ind3);
lon_inc = longs(ind3);
dep_inc = depths(ind3);
meanSR_inc = meanSR(ind3);
res_byage_inc = res_byage(ind3);
res_bydepth_inc = res_bydepth(ind3);

%Make plot with locations denoted as red stars
figure(23)
worldmap('World')
setm(gca, 'mapprojection', 'robinson')
geoshow('landareas.shp', 'FaceColor','[0.7 0.7 0.7]', 'EdgeColor', '[0.7 0.7 0.7]')
load coastlines
plotm(coastlat, coastlon, 'Color', 'k')
hold on
plotm(lat_inc, lon_inc, 'r*')
%% Plot Histograms

%Plot histograms of Mean SR, Depths, Resolution by Age, Resolution by Depth
figure(29)
subplot(4,1,1)
histogram(meanSR_inc, 0:2:90, 'FaceColor','k')
xlabel('Mean SR (cm/kyr)')
ylabel('Counts')
xlim([0 90])
title("Cores' Mean SR")

subplot(4,1,2)
histogram(dep_inc./1000, 0:0.25:6, 'FaceColor', 'k')
xlim([0 6])
xlabel("Depth (km)")
ylabel("Counts")
xticks(0:0.5:6)
title("Core Depths")

subplot(4,1,3)
histogram(res_byage_inc, 20, 'FaceColor','k')
xlabel('Sampling Frequency by age (kyr/date)')
ylabel('Counts')
title("Cores' Mean Resolution By Age")

subplot(4,1,4)
histogram(res_bydepth_inc, 20, 'FaceColor','k')
%xlim([0 150])
xlabel('Sampling Frequency by depth (cm/date)')
ylabel('Counts')
%xticks(0:25:150)
title("Cores' Mean Resolution By Depth")

%% Display some important information

%Input values from BIGMACS transition matrix
BIGMACSTM = [0.7215, 0.0940, 0.1846; 0.4328, 0.2687, 0.2985; 0.2670, 0.1041, 0.6290];

%Find out how much sediment the cores used constitute
lengthsed_core = nan(1,length(nSRcounts));
for i = 1:length(nSRcounts)
    if isempty(nSRcounts{i})
    else
    lengthsed_core(i) = sum(nSRcounts{i}(2,:));
    end
end

%Display total number fo 14C pairs used and total length of sediment used
disp("The total number of 14C pairs used is")
disp(sum(num14cpairs, 'omitmissing'))
disp("The total length of sediment used is")
disp(sum(sedimentlength, 'omitmissing'))
disp(sum(lengthsed_core, 'omitmissing'))

%results = struct("CoreNames", core_inc, "Latitude", lat_inc, "Longitude", lon_inc, "Depths", dep_inc, "NormSRDist", [SRvals_interp; SRvalsprob_norm], "TransitionMatrix", allTM, "HistogramBinCounts", bin, "HistogramBinEdges", binedges, "nSRcounts", nSRcounts, "MeanSR", meanSR_inc, "ResolutionByAge", res_byage, "ResolutionByDepth", res_bydepth, "Number14Cpairs", num14cpairs, "SedimentLength", sedimentlength);
