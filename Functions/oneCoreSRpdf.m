function [core_invSRvals, core_invSRprobs, meanSR, MSI_byage_mean, MSI_bydepth_mean, lengthsed_mean, numdatepairs_mean, ageModes, scenarios, label, numreversals_mean, scenario_meanSR] = oneCoreSRpdf(corename, LabIDs, incDepths, excLabIDs, excDepths, S, plotfigs)
%% Read in Radiocarbon Data
[age, depth_cm, error, label] = getDataWA(corename);

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
[uniqueLabels, uniqueIdx] = unique(label);
if length(uniqueIdx) ~= length(label);
    for ilabel = 1:length(label)
        label(ilabel) = "xFakeLabel" + num2str(ilabel);
    end
end

%Create scenarios to get around doubly-dated depths
[scenarios, duplicated_depths, chosenLabels] = scenariosDDD(depth_cm, label, corename);

%% For each scenario, run the calculations for calculating inverse sed rate
%Note, that if a scenario throws up a problematic age reversal, new
%scenarios are constructed to avoid this.
newscenarios = 1;
scenariosCFR = zeros(size(scenarios));

%Run scenarios deal with reversals until there are no more reversals in any
%scenarios
[scenarios2, scenarios2CFR, chosenLabels2, scenario_invSRvals, scenario_invSRprobs, scenario_meanSR, numreversals, numdatepairs, ageModes, lengthsed, newscenarios, MSI_byage, MSI_bydepth] = scenariosDealWithReversals(scenarios, scenariosCFR, chosenLabels, depth_cm, age, error, label, corename, duplicated_depths, S, plotfigs);
while newscenarios == 1
    %Check whether any scenarios are replicas of each other
    scenlengths = cellfun(@length, scenarios2); %Get lengths of each scenario
    Uscenlengths = unique(scenlengths);         %Find how many different lengths
    scenarios2keep = [];
    shifter = 0;
    for i = 1:length(Uscenlengths)
        scenLog = scenlengths == Uscenlengths(i);   %Find each scenario that is a given length
        sameLengthScens = scenarios2(scenLog);  %Get these scenarios in a cell array
        scenArray = strings(Uscenlengths(i), sum(scenLog));
        for j = 1:sum(scenLog)
            scenArray(:,j) = sameLengthScens{j}; %Convert them into a string array
        end
        [~, scenarios2keepI, ~] = unique(scenArray', 'rows'); %Find all the unique rows (akin to unique scenarios)
        if length(scenarios2keepI) ~= sum(scenLog)
            % disp("Getting rid of a duplicate scenario!")
        end
        scenarios2keep = [scenarios2keep; scenarios2keepI + shifter];
        shifter = shifter + sum(scenLog);
    end
    scenarios2 = scenarios2(scenarios2keep);
    scenarios2CFR = scenarios2CFR(scenarios2keep);
    chosenLabels2 = chosenLabels2(scenarios2keep);

    [scenarios2, scenarios2CFR, chosenLabels2, scenario_invSRvals, scenario_invSRprobs, scenario_meanSR, numreversals, numdatepairs, ageModes, lengthsed, newscenarios, MSI_byage, MSI_bydepth] = scenariosDealWithReversals(scenarios2, scenarios2CFR, chosenLabels2, depth_cm, age, error, label, corename, duplicated_depths, S,plotfigs);
end

scenarios = scenarios2;

%Calculate means of interested variables across all scenarios
meanSR = mean(scenario_meanSR, 'omitmissing');
numdatepairs_mean = mean(numdatepairs, 'omitmissing');
lengthsed_mean = mean(lengthsed, 'omitmissing');
numreversals_mean = mean(numreversals, 'omitmissing');
MSI_byage_mean = mean(MSI_byage, 'omitmissing');
MSI_bydepth_mean = mean(MSI_bydepth, 'omitmissing');

%Find all possible ageModes from all scenarios



%Combine the results from each scenario, to get a core-specific pdf of
%invSRvals
if ~isempty(scenarios)
[core_invSRvals, core_invSRprobs] = combinepdfs(scenario_invSRvals, scenario_invSRprobs, lengthsed);
else
    core_invSRvals = [];
    core_invSRprobs = [];
end

end