addpath('Functions')

%% Load Metadata of MSPF cores
%Check which cores have MSPF (monospecific planktonic foram) dates
data = readtable("COPYcorechoices_MSPF_highRes2.xlsx"); %read all metadata
dataMSPF = data(data.MSPF == 1,:);

%------ Index Good Cores
%Create index for all cores
numAllCores = length(dataMSPF.CoreName);
allcores = 1:numAllCores;
%Create index for a specific core
namedLog = contains(string(dataMSPF.CoreName), "MD01-2378");

% Get info from Excel Spreadsheet
chosenCoresLog = namedLog;
cores = table2array(dataMSPF(chosenCoresLog, "CoreName")); %take list of MSPF corenames
lats = table2array(dataMSPF(chosenCoresLog, "LatitudeDec"));
longs = table2array(dataMSPF(chosenCoresLog,"LongitudeDec"));
depths = table2array(dataMSPF(chosenCoresLog, "WaterDepthM"));
LabIDs = table2cell(dataMSPF(chosenCoresLog, "LabIDs")); %take list of LabIDs relating to MSPF dates of each core
incDepths = table2cell(dataMSPF(chosenCoresLog, "IncludeDepths")); % take list of depths (useful if no labels)
excLabIDs = table2cell(dataMSPF(chosenCoresLog, "excludeLabIDs")); %take list of manually removed dates for each core
excDepths = table2cell(dataMSPF(chosenCoresLog, "excludeDepth")); %take list of manually removed dates for each core (useful if no labels)
numCores = sum(chosenCoresLog);

% Plot all radiocarbon data
iPlot = 1;
corePlot(cores{iPlot}, LabIDs{iPlot}, incDepths{iPlot}, excLabIDs{iPlot}, excDepths{iPlot})

%Get reversal data, etc
i = 1;
[core_invSRvals{i}, core_invSRprobs{i}, meanSR(i), MSI_byage(i), MSI_bydepth(i), sedimentlength(i), num14cpairs(i), corescenarios{i}, newlabels{i}, numreversals(i)] = oneCoreSRpdf(cores{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, 0);

%Convert to SR
[SRvals_interp, SRvalsprob_norm] = invSRtoSR(core_invSRvals{i}, core_invSRprobs{i});

figure;
plot(SRvals_interp, SRvalsprob_norm)
xlim([0 6])
ylabel("Probability")
xlabel("Sed Rate Ratio")
legend("MSPF2 PDF")

calcTM = 0;
[cores_transnums(:,:,i), cores_CSE2x(:,:,i), nSRcounts{i}, agediffs{i}] = oneCoreTM(cores{i}, corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, calcTM);

plotSRandResHistograms(nSRcounts, agediffs, num14cpairs, 1, 101, 'k', cores{1})

%----- Get output metadata and summary figures for all cores used
%outputMetadataAndSummaryFigures(1, cores, lats, longs, depths, meanSR, MSI_byage, MSI_bydepth, nSRcounts, sedimentlength, num14cpairs)

% Go through process again but with restriction that age differences >
% 1000kyr
ii = 2;
[cores_transnums(:,:,ii), cores_CSE2x(:,:,ii), nSRcounts{ii}, agediffs{ii}] = oneCoreTMRestrict(cores{i}, corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, calcTM);
nSRcountsRestricted1000 = nSRcounts(2);
agediffsRestricted1000 = agediffs(2);

plotSRandResHistograms(nSRcountsRestricted1000, agediffsRestricted1000, num14cpairs, 1, 102, 'k', cores{1})
