function [core_invSRvals, core_invSRprobs, meanSR, MSI_byage_mean, MSI_bydepth_mean, lengthsed_mean, numdatepairs_mean, scenarios, label, numreversals_mean] = oneCoreSRpdf(corename, LabIDs, incDepths, excLabIDs, excDepths, plotfigs)
%% Read in Radiocarbon Data
[age, depth_cm, error, label] = getDataWA(corename);

%% Filtering
%Filter for MSPF dates, remove manually determined outliers, only keep
%dates between 1 and 42 14C ky BP, and only keep cores with 4 or more
%accepted dates

[age, depth_cm, error, label, emptybreak1, emptybreak2,~,~,~] = filtering(age, depth_cm, error, label, LabIDs, incDepths, excLabIDs, excDepths);

if emptybreak1 == 1 || emptybreak2 == 1
    core_invSRvals = [];
    core_invSRprobs = [];
    meanSR = NaN;
    MSI_byage_mean = NaN;
    MSI_bydepth_mean = NaN;
    numdatepairs_mean = NaN;
    lengthsed_mean = NaN;
    scenarios = [];
    numreversals_mean = NaN;
    return
end
    
%Plot radiocarbon ages that fit through filtering
if plotfigs == 1
    figure()
    errorbar(depth_cm, age, error, "vertical", 'o', "color", 'k')
    set(gca, 'YTickLabel',get(gca,'YTick'))
    xlabel("Depth (cm)")
    ylabel(["Radiocarbon Age","(14C kyr BP)"])
    hold on
end

%% Doubly Dated Depth Scenarios
%Create fake labIDs where there are none
if sum(contains(label, "NaN"))~=0
    for ilabel = find(contains(label, "NaN"))
        label(ilabel) = "FakeLabel" + num2str(ilabel);
    end
end

%Create scenarios to get around doubly-dated depths
[scenarios, duplicated_depths] = scenariosDDD(depth_cm, label, corename);

%% For each scenario, run the calculations for calculating inverse sed rate
%Note, that if a scenario throws up a problematic age reversal, new
%scenarios are constructed to avoid this.
newscenarios = 1;


%Run scenarios deal with reversals until there are no more reversals in any
%scenarios
[scenarios2, scenario_invSRvals, scenario_invSRprobs, scenario_meanSR, numreversals, numdatepairs, lengthsed, newscenarios, MSI_byage, MSI_bydepth] = scenariosDealWithReversals(scenarios, depth_cm, age, error, label, corename, duplicated_depths, plotfigs);
while newscenarios == 1
    [scenarios2, scenario_invSRvals, scenario_invSRprobs, scenario_meanSR, numreversals, numdatepairs, lengthsed, newscenarios, MSI_byage, MSI_bydepth] = scenariosDealWithReversals(scenarios2, depth_cm, age, error, label, corename, duplicated_depths, plotfigs);
end

scenarios = scenarios2;

%Calculate means of interested variables across all scenarios
meanSR = mean(scenario_meanSR, 'omitmissing');
numdatepairs_mean = mean(numdatepairs, 'omitmissing');
lengthsed_mean = mean(lengthsed, 'omitmissing');
numreversals_mean = mean(numreversals, 'omitmissing');
MSI_byage_mean = mean(MSI_byage, 'omitmissing');
MSI_bydepth_mean = mean(MSI_bydepth, 'omitmissing');


%Combine the results from each scenario, to get a core-specific pdf of
%invSRvals
[core_invSRvals, core_invSRprobs] = combinepdfs(scenario_invSRvals, scenario_invSRprobs, lengthsed);

end