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
% 

%% Add folder of necessary functions to path
addpath('Functions')

%% Load Metadata of MSPF cores
%Check which cores have MSPF (monospecific planktonic foram) dates
%rawdata     = readtable("COPYcorechoices_MSPF_highRes.xlsx"); %read all metadata
%sheet = "COPYcorechoices_MSPF_highRes.xlsx";
sheet = "COPYcore40Metadata.xlsx";
rawdata     = readtable(sheet);
if sheet == "COPYcore40Metadata.xlsx"
    rawdataMSPF = rawdata(rawdata.UseChoice == 1,:);
else
    rawdataMSPF = rawdata(rawdata.MSPF == 1,:);
end

%% Create settings structure
S.minNumberOfAges   = 4;      %Minimum number of ages a core must have (after filtering) to be used
S.DeltaRError       = 200;        % Error put on the Delta R (reservoir age correction)
S.reversalCriteria  = 0.75;  %What fraction of SRs between two ages must be negative to call it a reversal
S.weighting         = true;         %Whether or not to weight by depth (= false means no weighting)
S.normWithRunMean   = false;    %Use a common meanSR to use when calculating normalised SR (false) or use the SR from each individual run (true).

%% Choose cores to look at
%------ Index Good Cores
%Create index for all cores
numAllCores = length(rawdataMSPF.CoreName);
allcores    = 1:numAllCores;

%Create logical of cores without those that have been manually removed
%(goodLog)
reversalDenseCores = ["GeoB1711-4", "H214", "SO75_3_26KL", "KNR159-5-36GGC", "MD02-2550", "GIK17940-2"];
problemCores       = []; %reversals are mixed with duplicated depths, haven't figured out how to handle this yet
badLog             = ismember(string(rawdataMSPF.CoreName),[reversalDenseCores, problemCores]);
badIndexes         = allcores(badLog);
goodLog            = badLog == 0;
goodIndexes = allcores(~badLog);

%Create some other useful logical, to test all cores, only 1 core (defined by name), or a
%subset of cores (defined by index)
namedLog           = ismember(string(rawdataMSPF.CoreName), "A7");
allLog             = true(length(goodLog), 1);
subsetChooser      = logical(zeros(numAllCores, 1)); %#ok<LOGL>
subsetChooser(50:end)  = 1; 
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
if sheet ~= "COPYcore40Metadata.xlsx"

LabIDs      = table2cell(rawdataMSPF(chosenCoresLog, "LabIDs")); %take list of LabIDs relating to MSPF dates of each core
incDepths   = table2cell(rawdataMSPF(chosenCoresLog, "IncludeDepths")); % take list of depths (useful if no labels)
excLabIDs   = table2cell(rawdataMSPF(chosenCoresLog, "excludeLabIDs")); %take list of manually removed dates for each core
excDepths   = table2cell(rawdataMSPF(chosenCoresLog, "excludeDepth")); %take list of manually removed dates for each core (useful if no labels)


elseif sheet == "COPYcore40Metadata.xlsx"
    LabIDs = cell(numCores, 1);
    incDepths = cell(numCores, 1);
    excLabIDs = cell(numCores, 1);
    excDepths = cell(numCores, 1);
    %Set up all LabIDs to be "all"
    for i = 1:length(cores)
        LabIDs{i} = "all";
        incDepths{i} = "";
    end
    %Set up excLabIDs to be a combination of all reasons to exclude LabIDs
    reversalIDs = table2cell(rawdataMSPF(chosenCoresLog, "ReversalIDs"));
    AgeGapIDs = table2cell(rawdataMSPF(chosenCoresLog, "AgeGapIDs"));
    NonPlanktonicIDs = table2cell(rawdataMSPF(chosenCoresLog, "NonPlanktonicIDs"));
    MiscRemovalIDs = table2cell(rawdataMSPF(chosenCoresLog, "MiscRemovalIDs"));
    for i = 1:length(cores)
        IDvector = [string(reversalIDs{i}), string(AgeGapIDs{i}), string(NonPlanktonicIDs{i}), string(MiscRemovalIDs{i})];
        IDvector = IDvector(IDvector ~= "");
        IDstring = "";
        for j = 1:length(IDvector)
            if j == 1
                IDstring = IDvector(j);
            else
            IDstring = IDstring+ ", " + IDvector(j);
            end
        end
        %IDstring = [IDvector(1) + ", " + IDvector(2) + ", " + IDvector(3) + ", " + IDvector(4)];
        excLabIDs{i} = char(IDstring);
    end
    %Set up excDepths to be a combination of all reasons to exclude Depths
    ReversalDepths = table2cell(rawdataMSPF(chosenCoresLog, "ReversalDepths"));
    AgeGapDepths = table2cell(rawdataMSPF(chosenCoresLog, "AgeGapDepths"));
    for i = 1:length(cores)

        if isnan(ReversalDepths{i})
            ReversalDepths{i} = "";
        end
        DepthVector = [string(ReversalDepths{i}), string(AgeGapDepths{i})];
        DepthVector = DepthVector(DepthVector ~= "");
        DepthString = "";
        for j = 1:length(DepthVector)
            if j == 1
                DepthString = DepthVector(j);
            else
            DepthString = DepthString+ ", " + DepthVector(j);
            end
        end
        %DepthString = [DepthVector(1) + ", " + DepthVector(2) + ", " + DepthVector(3) + ", " + DepthVector(4)];
        excDepths{i} = char(DepthString);
    end

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
parfor i = 1:numCores
    disp("Working on " + cores{i})
    [core_invSRvals{i}, core_invSRprobs{i}, meanSR(i), MSI_byage(i),...
        MSI_bydepth(i), sedimentlength(i), num14cpairs(i), ageModes{i},...
        corescenarios{i}, newlabels{i}, numreversals(i), scenario_meanSR{i}]...
        = oneCoreSRpdf(cores{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, S, 0);
    disp("Done with " + cores{i})
end
disp("Created all scenarios")

% Find mean calibrated ages used in all scenarios of each core
% 
% for i = 1:numCores
%     scen = corescenarios{i};
%     usedIDs{i} = [];
%     for j = 1:length(scen)
%         scenIDs = scen{j};
%         newIDs = scenIDs(~contains(scenIDs, usedIDs));
%         usedIDs{i} = [usedIDs{i}; newIDs];
%     end
% 
% 
% end




%% Create a table that holds all useful metadata
agecoverage = sedimentlength./meanSR;
dataT = table(cores,lats, longs, depths, ocean, meanSR, ageModes, sedimentlength, agecoverage, num14cpairs, MSI_byage, MSI_bydepth, LabIDs, incDepths, excLabIDs, excDepths, core_invSRvals, core_invSRprobs, corescenarios);

%% ------ Define Subsets of interest
% high SR and low SR cores (separated by 8cm/kyr following Lin et al., 2014)
depth1000Log = depths > 1000;
lowSRCoresLog   = meanSR<= 8 & depth1000Log;
highSRCoresLog  = meanSR >8 & depth1000Log;
allCoresLog     = ~isnan(meanSR) & depth1000Log;

% find 10 highSR cores with most data
[~,highSRhighResCoresInd]   = maxk(num14cpairs.*highSRCoresLog, 10);
highSRhighResCoresLog       = unfind(highSRhighResCoresInd, numel(cores));

% %% plot figure showing all the ages being used
% figure;
% yTL = string(); %This will hold the core names in a format to be used for yticklabels
% %for i = 1:length(ageModes)
% indexes_hSR = find(highSRCoresLog);
% for i = 1:length(indexes_hSR)
%     %For each core, plot the mode of each calibrated radiocarbon date on a
%     %horizontal line, with a 1 pt vertical offset to compare cores
%     plot(ageModes{indexes_hSR(i)}./1000, i*ones(size(ageModes{indexes_hSR(i)})), 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'Color', 'r', 'LineWidth', 1)
%     hold on
%     yTL = [yTL, cores{indexes_hSR(i)}];
% end
% indexes_lSR = find(lowSRCoresLog);
% for i = 1:length(indexes_lSR)
%     %For each core, plot the mode of each calibrated radiocarbon date on a
%     %horizontal line, with a 1 pt vertical offset to compare cores
%     plot(ageModes{indexes_lSR(i)}./1000, length(indexes_hSR)+i*ones(size(ageModes{indexes_lSR(i)})), 'Marker', 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 8, 'Color', 'b', 'LineWidth', 1)
%     hold on
%     xlabel("Age (kyr)")
%     yTL = [yTL, cores{indexes_lSR(i)}];
% end
% yticks(1:length(yTL))
% yticklabels(yTL(2:end))
% ylim([0 length(yTL)+1])
% 
% %% plot figure showing all the ages being used
% figure;
% yTL = string(); %This will hold the core names in a format to be used for yticklabels
% %for i = 1:length(ageModes)
% indexes = find(depth1000Log);
% for i = 1:length(indexes)
%     %For each core, plot the mode of each calibrated radiocarbon date on a
%     %horizontal line, with a 1 pt vertical offset to compare cores
%     plot(ageModes{indexes(i)}./1000, i*ones(size(ageModes{indexes(i)})), 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'Color', 'k', 'LineWidth', 1)
%     hold on
%     xlabel("Age (kyr)")
%     yTL = [yTL, cores{indexes(i)}];
% end
% yticks(1:length(yTL))
% yticklabels(yTL(2:end))
% ylim([0 length(yTL)+1])

%% Simply make all ages, except those manually filtered, possible


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

parfor i =  1:numCores
    [nSRcounts{i}, agediffs{i}] = oneCoreTMRestrict(cores{i}, corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, scenario_meanSR{i}, S, 0);
end

% ------ Run through all chosen cores with random sampling approach,
% RESTRICTION ON MINIMUM AGE DIFFERENCE = 500

%Initiate variables
disp("b")
nSRcounts500   = cell(numCores, 1); % Holds all the nSR counts (which form histogram that makes nSR pdf)
agediffs500    = cell(numCores, 1); % Holds all the age differences for each nSR measurement (resolution pdf)

parfor i =  1:numCores
    [nSRcounts500{i}, agediffs500{i}] = oneCoreTMRestrict(cores{i}, corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, scenario_meanSR{i}, S, 500);
end

% ------ Run through all chosen cores with random sampling approach,
% RESTRICTION ON MINIMUM AGE DIFFERENCE = 1000

%Initiate variables
nSRcounts1000   = cell(numCores, 1); % Holds all the nSR counts (which form histogram that makes nSR pdf)
agediffs1000    = cell(numCores, 1); % Holds all the age differences for each nSR measurement (resolution pdf)

parfor i =  1:numCores
    [nSRcounts1000{i}, agediffs1000{i}] = oneCoreTMRestrict(cores{i}, corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, scenario_meanSR{i}, S, 1000);
end

% ------ Run through all chosen cores with random sampling approach,
% RESTRICTION ON MINIMUM AGE DIFFERENCE = 1000

%Initiate variables
nSRcounts1500   = cell(numCores, 1); % Holds all the nSR counts (which form histogram that makes nSR pdf)
agediffs1500    = cell(numCores, 1); % Holds all the age differences for each nSR measurement (resolution pdf)
% 
% parfor i =  1:numCores
%     [nSRcounts1500{i}, agediffs1500{i}] = oneCoreTMRestrict(cores{i}, corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, scenario_meanSR{i}, S, 1500);
% end

% ------ Run through all chosen cores with random sampling approach,
% RESTRICTION ON MINIMUM AGE DIFFERENCE = 1000

%Initiate variables
nSRcounts2000   = cell(numCores, 1); % Holds all the nSR counts (which form histogram that makes nSR pdf)
agediffs2000    = cell(numCores, 1); % Holds all the age differences for each nSR measurement (resolution pdf)

% parfor i =  1:numCores
%     [nSRcounts2000{i}, agediffs2000{i}] = oneCoreTMRestrict(cores{i}, corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, scenario_meanSR{i}, S, 2000);
% end

dataT = addvars(dataT, nSRcounts, nSRcounts500, nSRcounts1000, nSRcounts1500, nSRcounts2000); 

%save("Results/dataT_planktonicF17_Nov7.mat", "dataT", "S")

%% Calculate transition matrices for each setup
% TM0     = TMcalculation(nSRcounts);
% [~,~,TM1000]  = TMcalculation(nSRcounts1000, highSRCoresLog);

x = 0.01:0.01:15; plotSRandResHistograms(dataT.nSRcounts500, x, ones(size(dataT.nSRcounts500)), true, 3,1,2,0,"",true)