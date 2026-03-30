function [core_invSRvals, core_invSRprobs, meanSR, MSI_byage_mean, MSI_bydepth_mean, lengthsed_mean, numdatepairs_mean, ageModes, scenarios, label, numreversals_mean, scenario_meanSR] = oneCoreScenarios(corename, dataLoc, LabIDs, incDepths, excLabIDs, excDepths, S, plotfigs)
% oneCoreScenarios  Per-core preprocessing: filter dates, build age-depth
%               scenarios, handle reversals, and compute inverse SR PDF.
%
% This function is the main per-core computation called from calcData.m.
% For a single sediment core it:
%   1. Reads radiocarbon data from the specified source.
%   2. Filters dates (material type, manual exclusions, calibration limits).
%   3. Constructs age-depth scenarios to handle doubly-dated depths (DDDs).
%   4. Iteratively resolves age reversals within each scenario by splitting
%      or invalidating scenarios via scenariosDealWithReversals.
%   5. Computes summary statistics (mean SR, MSI, sediment length) across
%      all valid scenarios.
%   6. Optionally combines scenario inverse-SR PDFs into a single core PDF.
%
% The returned scenarios and ageModes are passed to oneCoreRSR for
% the RSR methods, and to nSRBchron for the Bchron methods.
%
% INPUTS
%   corename  - (string) Sediment core identifier (e.g. "RC13-228")
%   dataLoc   - (string) Radiocarbon data source: "WA" (World Atlas) or
%               "Lin2014"
%   LabIDs    - (string) Lab IDs to include (from extract3; passed to
%               filtering)
%   incDepths - (string) Depths to include (from extract3; passed to
%               filtering)
%   excLabIDs - (string) Lab IDs to exclude manually (passed to filtering)
%   excDepths - (string) Depths to exclude manually (passed to filtering)
%   S         - (struct) Settings struct. Key fields:
%                 .pdfMethod    (logical) Compute and combine inverse SR PDFs
%                 .minNumberOfAges  Minimum dates required
%                 (plus fields passed to filtering, scenariosDDD,
%                  scenariosDealWithReversals, multiMatcalQ)
%   plotfigs  - (logical) Whether to produce diagnostic plots
%
% OUTPUTS
%   core_invSRvals    - Inverse SR values from the combined core PDF
%                       (empty if S.pdfMethod is false or core rejected)
%   core_invSRprobs   - Probabilities corresponding to core_invSRvals
%   meanSR            - Mean SR across all valid scenarios (cm/kyr)
%   MSI_byage_mean    - Mean Spacing Index by age, averaged across scenarios
%   MSI_bydepth_mean  - Mean Spacing Index by depth, averaged across scenarios
%   lengthsed_mean    - Mean sediment length across scenarios (cm)
%   numdatepairs_mean - Mean number of consecutive date pairs per scenario
%   ageModes          - (cell array) Mode age for each date in each scenario
%   scenarios         - (cell array) Valid age-depth scenarios (each a
%                       string array of LabIDs); "invalid" for rejected ones
%   label             - (string vector) Lab IDs after filtering
%   numreversals_mean - Placeholder (currently returns 999; not yet
%                       implemented)
%   scenario_meanSR   - (numeric vector) Mean SR for each scenario (cm/kyr)
%
% See also: calcData, filtering, scenariosDDD, scenariosDealWithReversals,
%           oneCoreRSR, nSRBchron

%% Read in Radiocarbon Data
if dataLoc == "WA"
    [age, depth_cm, error, label] = getDataWA(corename, S);
elseif dataLoc == "Lin2014"
    [age, depth_cm, error, label] = getDatatxt(corename, S);
end

%% Filtering
%Filter for MSPF dates, remove manually determined outliers, only keep
%dates between 1 and 42 14C ky BP, and only keep cores with 4 or more
%accepted dates

[age, depth_cm, error, label, emptybreak1, emptybreak2,~,~,~] = filtering(age, depth_cm, error, label, LabIDs, incDepths, excLabIDs, excDepths, corename, S);

if emptybreak1 == 1 || emptybreak2 == 1
    core_invSRvals = [];
    core_invSRprobs = [];
    meanSR = NaN;
    MSI_byage_mean = NaN;
    MSI_bydepth_mean = NaN;
    numdatepairs_mean = NaN;
    ageModes = [];
    lengthsed_mean = NaN;
    scenarios = [];
    numreversals_mean = NaN;
    scenario_meanSR = [];
    return
end

%% Doubly Dated Depth Scenarios
%Create fake labIDs where there are none
if sum(contains(label, "NaN"))~=0
    for ilabel = find(contains(label, "NaN"))
        label(ilabel) = "FakeLabel" + num2str(ilabel);
    end
end

%Check to make sure there aren't duplicated labels in the core data. If
%there are, replace them with Fake Labels so they don't mess up tests for
%uniqueness later on
[~, uniqueIdx] = unique(label);
if length(uniqueIdx) ~= length(label)
    for ilabel = 1:length(label)
        label(ilabel) = "xFakeLabel" + num2str(ilabel);
    end
end

%Create scenarios to get around doubly-dated depths
[scenarios, duplicated_depths, chosenLabels] = scenariosDDD(depth_cm, label, corename);

%% Build scenarioStruct
% Bundle all per-scenario state into a single struct so it can be passed
% cleanly to scenariosDealWithReversals and removeDuplicateScenarios.
date_is = 1:length(age);
[ageprobAll, calAge] = multiMatcalQ(age, error, date_is, S);

numScenarios = length(scenarios);
SS.scenarios    = scenarios;
SS.CFR          = zeros(numScenarios, 1);   % 1 once a scenario is confirmed reversal-free
SS.chosenLabels = chosenLabels;
SS.invSRvals    = cell(numScenarios, 1);
SS.invSRprobs   = cell(numScenarios, 1);
SS.meanSR       = nan(numScenarios, 1);
SS.numdatepairs = nan(numScenarios, 1);
SS.ageModes     = cell(numScenarios, 1);
SS.lengthSed    = nan(numScenarios, 1);
SS.MSI_byage    = nan(numScenarios, 1);
SS.MSI_bydepth  = nan(numScenarios, 1);

IDpairs  = strings(0);
agediffV = {};

%% Iterative reversal resolution
% scenariosDealWithReversals processes scenarios one pass at a time. If it
% finds a reversal it modifies SS and returns newScenIndicator = 1, so we
% loop until every scenario has been confirmed reversal-free.
newScenIndicator = 1;
while newScenIndicator == 1
    SS = removeDuplicateScenarios(SS);

    [SS, newScenIndicator, IDpairs, agediffV] = ...
        scenariosDealWithReversals(SS, depth_cm, ageprobAll, calAge, ...
        label, corename, duplicated_depths, IDpairs, agediffV, S, plotfigs);

    if sum(size(SS.ageModes)) == 0
        a = 1; %#ok<NASGU>
    end
end

%% Unpack scenarioStruct for downstream calculations
scenarios       = SS.scenarios;
scenario_meanSR = SS.meanSR;
numdatepairs    = SS.numdatepairs;
lengthsed       = SS.lengthSed;
MSI_byage       = SS.MSI_byage;
MSI_bydepth     = SS.MSI_bydepth;
ageModes        = SS.ageModes;
scenario_invSRvals  = SS.invSRvals;
scenario_invSRprobs = SS.invSRprobs;

%% Summary statistics across all valid scenarios
meanSR            = mean(scenario_meanSR, 'omitmissing');
numdatepairs_mean = mean(numdatepairs, 'omitmissing');
lengthsed_mean    = mean(lengthsed, 'omitmissing');
numreversals_mean = 999; % placeholder — not yet implemented
MSI_byage_mean    = mean(MSI_byage, 'omitmissing');
MSI_bydepth_mean  = mean(MSI_bydepth, 'omitmissing');

%% Combine per-scenario PDFs into a single core PDF
if S.pdfMethod
    if ~isempty(scenarios)
        [core_invSRvals, core_invSRprobs] = combinepdfs(scenario_invSRvals, scenario_invSRprobs, lengthsed);
    else
        core_invSRvals = [];
        core_invSRprobs = [];
    end
else
    core_invSRvals = [];
    core_invSRprobs = [];
end

end
