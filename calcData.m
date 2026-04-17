%% calcData — Stage 1 of the nSRdist pipeline
%
% PURPOSE
%   Computes normalised sedimentation rate (nSR) histories for a chosen set
%   of marine sediment cores, using several independent estimation methods.
%   Results are saved to a .mat file in Results/ for use by fitData.m.
%
% PIPELINE OVERVIEW
%   calcData  →  (saves dataT, S, rawdataManual to Results/)
%   fitData   →  loads Results file, fits distributions, saves fit results
%
% NSR ESTIMATION METHODS
%   LinNSR   - Attempts an exact replication of Lin2014 method
%
%   BMedian  — Uses the median of the Bchron posterior age at each dated
%              depth to compute pairwise sedimentation rates; the mean SR
%              for normalisation is taken from the full Bchron model.
%              Conceptually analogous to the Lin et al. (2014) approach,
%              though not an exact replication.
%
%   BSamp    — Like BMedian but samples from the full Bchron posterior age
%              distribution (bchronProb) across S.numruns iterations.
%              Propagates age uncertainty more fully than BMedian.
%              Results stored in bchronProb cell array.
%
%   RSR0 / RSR500 / RSR1000 / RSR1500
%            — Random Sampling from calibrated Radiocarbon (RSR) age PDFs,
%              with a minimum age-difference restriction (dt_min = 0, 500,
%              1000, or 1500 yr) between consecutive dated depths. Uses
%              oneCoreRSR.m; the restriction filters out pairs too close in
%              age to give reliable sedimentation rate estimates.
%
% DOUBLY-DATED DEPTHS AND REVERSALS
%   Before RSR nSR estimation, oneCoreScenarios.m is called for each core.
%   It builds "scenarios" — all valid orderings of radiocarbon dates that
%   avoid age reversals — by handling:
%     - Doubly-dated depths (DDDs): multiple dates at the same depth
%     - Age reversals: probabilistic reversals are resolved by constructing
%       alternative scenario pathways through the dated sequence
%   The resulting scenarios are shared by all RSRx methods and inform the
%   per-scenario mean SR used for normalisation.
%
% INPUT DATA
%   Data are loaded from an Excel spreadsheet (S.sheet) which records, for
%   each core:
%     - Core name, location, water depth, ocean basin
%     - Which radiocarbon dates to include (planktonic foraminifera only)
%     - Which dates to exclude (benthic foram, outliers, large gaps)
%     - Flags linking rows to the Mulitza 2022 World Atlas (useLin/usePF)
%   Raw radiocarbon measurement data are read from the Mulitza 2022 dataset
%   (S.WApath) or from the Lin2014 supplementary database.
%
% OUTPUT (saved to Results/S.dataOutputFile.mat)
%   dataT        — Table with one row per core containing:
%                    cores, lats, longs, depths, ocean, meanSR, ageModes,
%                    sedimentlength, agecoverage, num14cpairs,
%                    MSI_byage, MSI_bydepth,
%                    LabIDs, incDepths, excLabIDs, excDepths, dataLoc,
%                    core_invSRvals, core_invSRprobs,
%                    corescenarios, scenario_meanSR,
%                    RSR0, RSR500, RSR1000, RSR1500,
%                    BMode, BMedian, BProb, LinNSR
%   S            — Settings structure used for this run
%   rawdataManual — Subset of the raw metadata table used
%
% SETTINGS (S structure — set in the "Create settings structure" section)
%   See inline comments on each field below.
%
% DEPENDENCIES
%   Functions/      — all helper functions (added to path at startup)
%   MatCal          — matcal / matcalq for radiocarbon calibration
%   R + Bchron      — called via Rscript for Bchron age-depth modelling

%% Add folder of necessary functions to path
addpath('Functions')

%% Create settings structure
% All run parameters are stored in S and saved alongside the results, so
% any output file is self-documenting.

% --- Paths ---
scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    error('Please run calcData.m as a script file, not pasted into the Command Window.');
end
S.sandboxPath = fileparts(scriptPath); % Root directory of the repository
S.RscriptPath = "/usr/local/bin/Rscript";  % Path to Rscript; find yours with `which Rscript` in Terminal
S.WApath       = "WA_Foraminiferal_Isotopes_2022";         % Path to the Mulitza 2022 World Atlas dataset
S.sheet        = "DataSheets/datasheet.xlsx"; % Metadata spreadsheet controlling which cores and dates are used

% --- Output file ---
stringID = "Run1";
S.dataOutputFile = "dataT_" + stringID; % Saved to Results/S.dataOutputFile.mat

% --- Core selection filters ---
S.minimumCoreDepth = 1000; % Minimum water depth (m); excludes shallow-water cores
S.maxAtlanticLatN   = 40;  % Maximum northern latitude for Atlantic cores (degrees)
S.maxOtherLatN      = 40;  % Maximum northern latitude for non-Atlantic cores
S.maxAtlanticLatS   = 40;  % Maximum southern latitude for Atlantic cores (degrees)
S.maxOtherLatS      = 40;  % Maximum southern latitude for non-Atlantic cores

% --- Data source and pre-processing ---
S.replicateLin2014  = false; % If true, overrides other settings to match Lin2014 exactly (see block below)
S.useLin            = true;  % Include cores from the Lin2014 database
S.modifyLin2014Data = true;  % If true, applies additional reversals / gap removal on Lin2014 dates
S.usePF             = true;  % Include cores added from the Mulitza 2022 dataset (planktonic foram dates only)
S.DeltaRError       = 200;   % 1-sigma error (yr) on the marine reservoir age correction (DeltaR); propagated through calibration
S.c14AgeLim         = [0 50];% Accepted radiocarbon age range (kyr); dates outside this window are excluded
S.normWithRunAve    = true;  % If true, normalises each RSR run's SRs by that run's own mean SR;
                             %   if false, uses a single shared mean SR across all runs for a core
S.numruns           = 1000;  % Number of Monte Carlo iterations for BSamp and RSRx methods

% --- Bchron settings ---
S.useBchron              = true;    % If false, skips Bchron approaches entirely
S.BchronFilter           = true;    % If true, removes obvious outliers and large age gaps before running Bchron
S.BchronFolderName       = "All1_RLGtrue_DS0p05_8Apr26"; % Subfolder under BchronInputs/ from which to read pre-computed Bchron outputs
S.BchronOutlier          = 0.05;    % OutlierProbs value passed to Bchron (prior probability of each date being an outlier)
S.BchronReversalCriteria = 0.75;    % Fraction of Bchron MCMC samples that must reject an age for it to be flagged as a reversal
S.BchronCalCurve         = "Marine20"; % Calibration curve to use in Bchron ("Marine20" or "Marine09" for Lin2014 replication)
S.BchronDepthSpacing     = 0.05;    % Depth interval (cm) at which Bchron reports age estimates along the core
S.BchronReDo             = false;   % If false, re-uses existing Bchron outputs if available; if true, reruns Bchron for all cores

% --- RSRx method settings ---
S.minNumberOfAges   = 4;     % Minimum number of accepted dates a core must have to be processed
S.reversalCriteria  = 0.75;  % Fraction of sampled SR estimates between two adjacent dates that must be
                             %   negative to classify the pair as a reversal and trigger scenario branching
S.removeLargeGaps   = true;  % If true, manually excludes date pairs with very large age gaps before scenario construction
S.pdfMinVal         = 1e-6;  % Probability density cutoff below which calibrated age PDF tails are truncated,
                             %   reducing vector length and accelerating calculations
S.pdfMethod         = false; % If true, runs the full analytical PDF method for all date pairs;
                             %   if false, uses PDF only for DDDs and reversals, sampling otherwise
S.useModes          = false; % If true, uses the mode of each calibrated age PDF as a point estimate
                             %   rather than sampling (diagnostic / fast-approximate mode)

% --- Lin2014 approach
S.Lin2014Method     = false; % Whether to calculate NSRs with Lin2014 method  

% --- Other settings
S.matcalFast = true;         % Whether to use my modified version of matcal. Removes need to add MatCal to matlab (copy stored within repo). Speeds up scenario construction for RSRx methods by ~10x.
S.plotAgeModes = false;      % Whether to plot all age modes after creation of scenarios
S.plotAgeDepthModels = false;% Whether to plot calibrated ages against depth for each core
%% Do I want to use Lin2014 set up?
if S.replicateLin2014 == 1
    S.modifyLin2014Data = false;
    S.DeltaRError = 0;
    S.c14AgeLim = [0 55];
    S.useBchron = true;
    S.BchronCalCurve = "Marine09";
    S.BchronReDo = false;

    %Set up core choice settings (first two explicitly stated, others
    %implicitly taken from cores they used)
    S.minimumCoreDepth = 1000;%Minimum core depth (mbsf) in m
    S.maxAtlanticLatN   = 40; %Maximum north latitude of Atlantic Cores (explicitly stated in paper)
    S.maxOtherLatN      = 52; % (not explicitly stated in paper, max LatN is 51ish)
    S.maxAtlanticLatS   = 50; % (not explicitly stated in paper, max LatS is -46ish)
    S.maxOtherLatS      = 50; % (not explicitly stated in paper, max LatS is -46ish)
end
%% Load metadata from the spreadsheet
% The spreadsheet (S.sheet) has one row per core and records:
%   CoreName, LatitudeDec, LongitudeDec, WaterDepthM, Basin
%   Lin2014    — 1 if the core was used by Lin et al. (2014)
%   Lin2014Keep — 1 if the core should be kept for replication even if not in WAuse
%   WAuse      — 1 if the core should be used from the Mulitza 2022 dataset
% For each core, LabIDs of included dates (planktonic foram only) and
% excluded dates (benthic foram, outliers, large gaps) are also stored.
rawdata     = readtable(S.sheet);

% Select relevant information based on what subsets of data has been chosen
if S.useLin && ~S.usePF
    rawdataManual = rawdata(rawdata.Lin2014 == 1, :);
elseif ~S.useLin && S.usePF
    rawdataManual = rawdata(rawdata.WAuse == 1 & isnan(rawdata.Lin2014), :);
elseif S.useLin && S.usePF
    rawdataManual = rawdata(rawdata.WAuse == 1 | rawdata.Lin2014Keep == 1,:);
end

%% Create logical to choose all good cores
%Number of All Cores
numAllCores = length(rawdataManual.CoreName);

%Name of any problem cores (useful if wanting to exclude a single core)
problemCores = "KNR197-3-36GGC";
badLog       = ismember(string(rawdataManual.CoreName),problemCores);

%Restrict which cores to analyse based on metadata, such as depth or
%latitude
AtlanticLatLog = contains(string(rawdataManual.Basin), "Atlantic") & rawdataManual.LatitudeDec <= S.maxAtlanticLatN & rawdataManual.LatitudeDec >= -S.maxAtlanticLatS;
OtherBasinLatLog = ~contains(string(rawdataManual.Basin), "Atlantic") & rawdataManual.LatitudeDec <= S.maxOtherLatN & rawdataManual.LatitudeDec >= -S.maxOtherLatS;
latitudeRestrictionLog = AtlanticLatLog | OtherBasinLatLog;
depthRestrictionLog    = rawdataManual.WaterDepthM >= S.minimumCoreDepth;

%Therefore create logical of all good cores
restrictions       = latitudeRestrictionLog & depthRestrictionLog;
goodLog            = badLog == 0 & restrictions == 1;

%% Create some other useful logicals
%Test a core based on it's name
 namedLog = ismember(string(rawdataManual.CoreName), "M35003-4");   

%Test a subset of cores
subsetChooser = false(numAllCores,1);
subsetChooser(101:119) = true;

%Test all cores
% subsetChooser = true(numAllCores, 1);

%% Get material data from Excel
% Get relevant metadata from Excel Spreadsheet into useful variables
chosenCoresLog = goodLog;
numCores    = sum(chosenCoresLog);
cores       = table2array(rawdataManual(chosenCoresLog, "CoreName")); %take list of MSPF corenames
lats        = table2array(rawdataManual(chosenCoresLog, "LatitudeDec"));
longs       = table2array(rawdataManual(chosenCoresLog, "LongitudeDec"));
depths      = table2array(rawdataManual(chosenCoresLog, "WaterDepthM"));
ocean       = table2array(rawdataManual(chosenCoresLog, "Basin"));
%distance2coast = table2array(rawdataManual(chosenCoresLog, "DistanceToCoast"));

[LabIDs, incDepths, excLabIDs, excDepths, dataLoc] = extract3(rawdataManual, chosenCoresLog, S);

rawdataUse = rawdataManual(chosenCoresLog, :);
%% Plot calibrated radiocarbon dates against depth
if S.plotAgeDepthModels
 for iPlot = (1:length(cores))
  corePlotCal(cores{iPlot}, LabIDs{iPlot}, incDepths{iPlot}, excLabIDs{iPlot}, excDepths{iPlot}, dataLoc(iPlot), S)
 end
end

%% Compute nSR histories using Bchron age-depth models (BMode, BMedian, BSamp)
% nSRBchron runs Bchron for each core (or loads existing results) and
% returns nSR histories in 3-row matrix format for each method:
%   BMode   — nSR from the modal Bchron age at each dated depth
%   BMedian — nSR from the median Bchron age at each dated depth
%   BProb   — nSRs sampled from the Bchron posterior (S.numruns samples)
%   LinNSR    — nSR history matching the Lin2014 mean SR normalisation
% Skipped entirely if S.useBchron is false.
BMode   = cell(numCores, 1);
BMedian = cell(numCores, 1);
BSamp   = cell(numCores, 1);
LinNSR    = cell(numCores, 1);
bchronMeanSRs= NaN(numCores, 1);

if S.useBchron
    for i = 1:numCores
        [BMode{i}, BMedian{i}, BSamp{i}, LinNSR{i},  bchronMeanSRs(i)] = nSRBchron(cores{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, dataLoc(i), S);
    end
end
%% Build scenarios and compute per-core summary statistics (oneCoreScenarios)
% oneCoreScenarios is called for every core before the RSRx sampling begins.
% It serves two purposes:
%   1. Scenario construction — resolves doubly-dated depths (DDDs) and age
%      reversals by enumerating all valid date orderings ("scenarios"). Each
%      scenario is a specific assignment of ages to depths that produces no
%      reversals. These scenarios are later used by oneCoreRSR.
%
%   2. Summary statistics — computes the mean SR, mean sampling interval
%      (MSI by age and by depth), total sediment length, number of 14C date
%      pairs, per-depth age modes, and per-scenario mean SRs used for
%      normalisation in the RSRx methods.
%
% Also computes an inverse-SR PDF (core_invSRvals / core_invSRprobs) via
% the analytical PDF method for diagnostic purposes (stored in dataT).

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

%Calculate SR distribution for each core, as well as meanSR and other
%useful information
tic
for i = 1:numCores
     disp("Working on " + cores{i})
    [core_invSRvals{i}, core_invSRprobs{i}, meanSR(i), MSI_byage(i),...
        MSI_bydepth(i), sedimentlength(i), num14cpairs(i), ageModes{i},...
        corescenarios{i}, newlabels{i}, numreversals(i), scenario_meanSR{i}]...
        = oneCoreScenarios(cores{i}, dataLoc(i), LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, S, 0);
end
toc
disp("Created all scenarios")

%% Create a table that holds all useful metadata
agecoverage = sedimentlength./meanSR;
dataT = table(cores,lats, longs, depths, ocean, meanSR, ageModes,...
    sedimentlength, agecoverage, num14cpairs, MSI_byage, MSI_bydepth,...
    LabIDs, incDepths, excLabIDs, excDepths, dataLoc, core_invSRvals,...
    core_invSRprobs, corescenarios, scenario_meanSR);

%% Plot figure showing all the ages being used
if S.plotAgeModes
    plotAgeModes(chosenCoresLog, chosenCoresLog, ageModes, cores)
end

%% Compute nSR histories by Random Sampling from calibrated Radiocarbon PDFs (RSRx)
% oneCoreRSR samples S.numruns times from the calibrated radiocarbon age
% PDFs of each core and computes pairwise sedimentation rates. The output
% is a 3-row nSR matrix (see README "nSR matrix format") holding all sampled
% nSR values and their associated depth/age differences.
%
% Four variants are run, differing only in the minimum age difference (dt_min)
% required between consecutive dated depths:
%   RSR0    — no restriction (any dt accepted)
%   RSR500  — pairs with dt < 500 yr excluded
%   RSR1000 — pairs with dt < 1000 yr excluded
%   RSR1500 — pairs with dt < 1500 yr excluded
%
% The restriction reduces bias from closely-spaced dates whose calibrated
% age PDFs overlap substantially. RSR500 and above use parfor for speed.
%
% For each variant, nSRcounts holds the nSR matrix cells and agediffs holds
% the corresponding age-difference histograms (used to characterise the
% resolution distribution of the dataset).

%------ RSR0: no minimum age-difference restriction

%Initiate variables
disp("Starting RSR0")
RSR0       = cell(numCores, 1); % Holds RSR0 data

for i =  1:numCores
    disp(cores{i})
    [RSR0{i}, ~] = oneCoreRSR(cores{i}, dataLoc(i), corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, scenario_meanSR{i}, ageModes{i}, S, 0);
end

%------ RSR500: minimum age-difference restriction = 500 yr

%Initiate variables
disp("Starting RSR500")
RSR500   = cell(numCores, 1); % Holds RSR500

parfor i =  1:numCores
    [RSR500{i}, ~] = oneCoreRSR(cores{i}, dataLoc(i), corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, scenario_meanSR{i}, ageModes{i}, S, 500);
end

disp("Starting RSR1000")
%------ RSR1000: minimum age-difference restriction = 1000 yr

%Initiate variables
RSR1000   = cell(numCores, 1); % Holds RSR1000

parfor i =  1:numCores
    [RSR1000{i}, ~] = oneCoreRSR(cores{i}, dataLoc(i), corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, scenario_meanSR{i}, ageModes{i}, S, 1000);
end
disp("Starting RSR1500")
%------ RSR1500: minimum age-difference restriction = 1500 yr

%Initiate variables
RSR1500   = cell(numCores, 1); % Holds all the nSR counts (which form histogram that makes nSR pdf)

parfor i =  1:numCores
    [RSR1500{i},~] = oneCoreRSR(cores{i}, dataLoc(i), corescenarios{i}, LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, scenario_meanSR{i}, ageModes{i}, S, 1500);
end
disp("All RSRx samples finished")

%% Assemble final table and save
% Append all nSR matrix cell arrays and Bchron outputs to dataT, then save
% everything (dataT, S, rawdataManual) to Results/S.dataOutputFile.mat.
% This file is the input to fitData.m.
dataT = addvars(dataT, RSR0, RSR500, RSR1000, RSR1500, BMode, BMedian, BSamp, LinNSR);
save(fullfile("Results", S.dataOutputFile), "dataT", "S", "rawdataManual")

