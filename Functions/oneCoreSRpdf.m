function [core_invSRvals, core_invSRprobs, meanSR, MSI_byage_mean, MSI_bydepth_mean, lengthsed_mean, numdatepairs_mean, ageModes, scenarios, label, numreversals_mean, scenario_meanSR] = oneCoreSRpdf(corename, dataLoc, LabIDs, incDepths, excLabIDs, excDepths, S, plotfigs)
%% Read in Radiocarbon Data
if dataLoc == "WA"
    [age, depth_cm, error, label] = getDataWA(corename);
elseif dataLoc == "Lin2014"
    [age, depth_cm, error, label] = getDatatxt(corename);
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
    
% %Plot radiocarbon ages that fit through filtering
% if plotfigs == 1
%     figure()
%     errorbar(depth_cm, age, error, "vertical", 'o', "color", 'k')
%     set(gca, 'YTickLabel',get(gca,'YTick'))
%     xlabel("Depth (cm)")
%     ylabel(["Radiocarbon Age","(14C kyr BP)"])
%     hold on
% end

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

%% For each scenario, run the calculations for calculating inverse sed rate
%Note, that if a scenario throws up a problematic age reversal, new
%scenarios are constructed to avoid this.

date_is = 1:length(age);
[ageprobAll, calAge] = multiMatcalQ(age, error, date_is, S);

%Run scenarios deal with reversals until there are no more reversals in any
%scenarios
scenariosCFR = zeros(size(scenarios)); %Initialise vector to signify whether a given scenario has been confirmed to have no reversals (1 if no reversals, 0 if not yet) 
IDpairs = strings(0);
agediffV = {};
numScenarios = length(scenarios);
scenario_invSRvals  = cell(numScenarios, 1);
scenario_invSRprobs = cell(numScenarios, 1);
scenario_meanSR     = nan(numScenarios, 1);
numdatepairs        = nan(numScenarios, 1);
ageModes            = cell(numScenarios, 1);
lengthsed           = nan(numScenarios, 1);
MSI_bydepth         = nan(numScenarios, 1);
MSI_byage           = nan(numScenarios, 1);

[scenarios, scenariosCFR, chosenLabels, scenario_invSRvals,...
    scenario_invSRprobs, scenario_meanSR, numdatepairs,...
    ageModes, lengthsed, newscenarios, MSI_byage, MSI_bydepth, IDpairs, agediffV]...
    = scenariosDealWithReversals(scenarios, scenariosCFR, chosenLabels,... %%%%%%
    scenario_invSRvals, scenario_invSRprobs, scenario_meanSR, ...
    numdatepairs,ageModes, lengthsed, MSI_byage, MSI_bydepth, ...
    depth_cm, ageprobAll, calAge, label, corename, duplicated_depths,...
    IDpairs, agediffV, S, plotfigs);

% if newscenarios == 1
%     scenariosCFR(:) = 0; %This ensures all scenarios get ageMode, lengthsed, etc information
% end

while newscenarios == 1
    %Check whether any scenarios are replicas of each other
    scenlengths = cellfun(@length, scenarios); %Get lengths of each scenario
    Uscenlengths = unique(scenlengths);         %Find how many different lengths
    scenarios2keep = [];
    shifter = 0;
    for i = 1:length(Uscenlengths)
        scenLog = scenlengths == Uscenlengths(i);                           %Find each scenario that is a given length
        sameLengthScens = scenarios(scenLog);                              %Get these scenarios in a cell array
        scenArray = strings(Uscenlengths(i), sum(scenLog));
        for j = 1:sum(scenLog)
            scenArray(:,j) = sameLengthScens{j};                            %Convert them into a string array
        end
        [~, scenarios2keepI, ~] = unique(scenArray', 'rows');               %Find all the unique rows (akin to unique scenarios)
        scenarios2keepIs = sort(scenarios2keepI);
        scenarios2keep = [scenarios2keep; scenarios2keepIs + shifter];       %#ok<AGROW>
        shifter = shifter + sum(scenLog);
    end
    scenarios = scenarios(scenarios2keep);
    scenariosCFR = scenariosCFR(scenarios2keep);
    chosenLabels = chosenLabels(scenarios2keep);
    scenario_invSRvals = scenario_invSRvals(scenarios2keep);
    scenario_invSRprobs = scenario_invSRprobs(scenarios2keep);
    scenario_meanSR = scenario_meanSR(scenarios2keep);
    numdatepairs = numdatepairs(scenarios2keep);
    ageModes =  ageModes(scenarios2keep);
    lengthsed = lengthsed(scenarios2keep);
    MSI_byage = MSI_byage(scenarios2keep);
    MSI_bydepth = MSI_bydepth(scenarios2keep);

    % % disp(num2str(sum(scenariosCFR)) + " scenarios confirmed with no reversals")
      %disp("size of scenarios vector " + num2str(length(scenarios)))

    [scenarios, scenariosCFR, chosenLabels, scenario_invSRvals,...
        scenario_invSRprobs, scenario_meanSR, numdatepairs,...
        ageModes, lengthsed, newscenarios, MSI_byage, MSI_bydepth, IDpairs, agediffV]...
     = scenariosDealWithReversals(scenarios, scenariosCFR, chosenLabels,... %%%%%%
    scenario_invSRvals, scenario_invSRprobs, scenario_meanSR, ...
    numdatepairs,ageModes, lengthsed, MSI_byage, MSI_bydepth, ...
    depth_cm, ageprobAll, calAge, label, corename, duplicated_depths,...
    IDpairs, agediffV, S, plotfigs);

    if sum(size(ageModes)) == 0;
        a = 1;
    end
end

%Calculate means of interested variables across all scenarios
meanSR = mean(scenario_meanSR, 'omitmissing');
numdatepairs_mean = mean(numdatepairs, 'omitmissing');
lengthsed_mean = mean(lengthsed, 'omitmissing');
numreversals_mean = 999; %mean(numreversals, 'omitmissing'); %Making this 999 for ease right now
MSI_byage_mean = mean(MSI_byage, 'omitmissing');
MSI_bydepth_mean = mean(MSI_bydepth, 'omitmissing');

%Find all possible ageModes from all scenarios

%Combine the results from each scenario, to get a core-specific pdf of
%invSRvals
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