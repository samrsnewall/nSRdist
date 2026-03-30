%Set up necessary paths
addpath('../../Functions')

%Set up necessary parameters
S.DeltaRError       = 200;      %Error put on the Delta R (reservoir age correction)
S.c14AgeLim         = [0 50];   %Cutoffs for radiocarbon ages, in kyr
S.minNumberOfAges   = 4;        %Minimum number of ages a core must have (after filtering) to be used
S.reversalCriteria  = 0.75;     %What fraction of SRs between two ages must be negative to call it a reversal
S.removeLargeGaps   = true;     %Whether to manually remove large age gaps or leave them in
S.pdfMinVal         = 1e-6;     %Cutoff to reduce size of radiocarbon pdf vector to accelerate calculations
S.pdfMethod         = false;    %Whether to complete all of the pdf method or simply use it for multiply dated depths and reversals
S.useModes          = false;    %Calculate nSR distribution with the mode of each radiocarbon distribution, not sampling

%Set up example radiocarbon set
label = ["ID1", "ID2", "ID3", "ID4", "ID5", "ID6", "ID7", "ID8"]';
depth_cm = [10, 15, 20, 20, 23, 25, 28, 31, 35]';
age = [2, 3.2, 4, 4.3, 3, 5.1, 5.9, 6.1, 6.5]';
error = [0.050, 0.050, 0.050, 0.060, 0.060, 0.060, 0.080, 0.080]';

figure;
subplot(2,1,1)
errorbar(depth_cm, age, error, 'vertical', 'LineStyle', 'none', 'Color', 'k')

%Set up info that would come from data sheet
LabIDs = "all";
incDepths = [];
excLabIDs = [];
excDepths = [];
corename = "FakeCore1";

%Perform filtering
[age, depth_cm, error, label, emptybreak1, emptybreak2,~,~,~] = filtering(age, depth_cm, error, label, LabIDs, incDepths, excLabIDs, excDepths, corename, S);

%Create scenarios to get around doubly-dated depths
[scenarios, duplicated_depths, chosenLabels] = scenariosDDD(depth_cm, label, corename);

%Calibrate ages
date_is = 1:length(age);
[ageprobAll, calAge] = multiMatcalQ(age, error, date_is, S);

%Build scenarioStruct and run iterative reversal resolution
numScenarios = length(scenarios);
SS.scenarios    = scenarios;
SS.CFR          = zeros(numScenarios, 1);
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
plotfigs = false;

newScenIndicator = 1;
while newScenIndicator == 1
    SS = removeDuplicateScenarios(SS);
    [SS, newScenIndicator, IDpairs, agediffV] = ...
        scenariosDealWithReversals(SS, depth_cm, ageprobAll, calAge, ...
        label, corename, duplicated_depths, IDpairs, agediffV, S, plotfigs);
end

%Unpack scenarioStruct
scenarios = SS.scenarios;
ageModes  = SS.ageModes;

subplot(2,1,2)
S1_logi = contains(label, scenarios{1});
hold on
errorbar(depth_cm, age, error, 'vertical', 'LineStyle', 'none', 'Color', 'k')

errorbar(depth_cm(S1_logi), age(S1_logi), error(S1_logi), 'vertical')

S2_logi = contains(label, scenarios{2});
%plot(depth_cm(S2_logi), age(S2_logi), '.r')
errorbar(depth_cm(S2_logi), age(S2_logi), error(S2_logi), 'vertical')

%So that sets up the scenarios from some data. Now I need to do the RSR
%part

