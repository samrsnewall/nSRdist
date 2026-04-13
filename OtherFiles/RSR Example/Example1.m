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
label = ["ID1", "ID2", "ID3", "ID4", "ID5", "ID6", "ID7", "ID8", "ID9"]';
depth_cm = [10, 15, 20, 20, 23, 25, 28, 31, 35]';
age = [2, 3.2, 4, 4.3, 3, 5.1, 5.9, 6.1, 6.5]';
error = [0.050, 0.050, 0.050, 0.060, 0.060, 0.060, 0.080, 0.080, 0.080]';

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

%% ------- Check which scenarios to analyse
%Set up empty vectors
validScenariosBool = ones(length(scenarios),1);

%Check which scenarios are invalid (due to reversals etc).
for i_sce = 1:length(scenarios)
    if scenarios{i_sce} == "invalid"
        validScenariosBool(i_sce) = 0;
    end
end

%% ------- Loop through scenarios using random sampling SR calculations
minAgeDiff = 1000;
S.normWithRunAve = true;

%Perform Calibrations
[ageprobAll, calAge] = multiMatcalQ(age, error, 1:length(age), S);

% Randomly choose which scenario to use for each run
numruns = 2;
rndScen = [1, 2];

%

%Create numruns many possible sedimentation histories, by random sampling
for ix = 1:numruns
    %Find which dates are in the random scenario
    date_bool = ismember(string(label), scenarios{rndScen(ix)});
    date_is = find(date_bool == 1);

    %Get ages' probability vectors
    ageprob = ageprobAll(:,date_is);

    %get depths in this scenario
    depths_scenario = depth_cm(date_is);

    %Assign single dates from the probability vectors either...
    %%%% Using modes
    if S.useModes
        run_age = ageModes{rndScen(ix)};
        run_age = run_age';
    else

        %%%% Using random sampling

        %Initialise
        run_age = NaN(1,length(date_is));

        %Run through each date and randomly choose a year from it's
        %probability vector, ensuring that each age is older than the age
        %directly shallower than it
        i = 1;
        while i <= length(date_is)
            %Set up age probabilities so that every age is older than the
            %age of sample shallower than it
            if i == 1
                ind1        = ageprob(:,i) >= S.pdfMinVal;
                poss_ages   = calAge(ind1);
                age_probs   = ageprob(ind1,i);
                run_age(i)  = randsample(poss_ages, 1, true, age_probs);
            else
                min_age = max(run_age(1:i-1));
                possibleAgesLog = calAge>min_age;
                poss_ages = calAge(possibleAgesLog);
                age_probs = ageprob(possibleAgesLog,i);
                age_probs = age_probs./sum(age_probs);
                ind1 = age_probs >= S.pdfMinVal;
                poss_ages = poss_ages(ind1);
                age_probs = age_probs(ind1);
                if length(poss_ages) <= 1
                    disp("Random sampling has reached the end of the ages vector (from Matcal) for core "  + corename + ": - next random sample will fail")
                    disp("Taking precaution - forcing assumption that no radiocarbon ages are actually older than 55000 cal yr")
                    %resetting this random sampling run
                    i = 1;
                    run_age = NaN(1,length(date_is));
                    continue
                end

                run_age(i) = randsample(poss_ages, 1, true, age_probs);
                %%% RESTRICTION CODE
                %If next age is not more than x years greater than
                %previous accepted age, reject it and move to sample from next
                %radiocarbon date (x stored in input minAgeDiff)
                if run_age(i) - max(run_age(1:i-1)) < minAgeDiff
                    run_age(i) = NaN;
                end
            end
            i = i+1;
        end
    end

    % Find out which ages were used and which weren't
    datesUsed = ~isnan(run_age);
    depthOfUsed = depths_scenario(datesUsed);
    runAgesOfUsed = run_age(datesUsed); %y
    %datesUsedStore(runN,:) = datesUsed;

    % Calculate average SR for that potential run
    aveSR_run = (depthOfUsed(end)-depthOfUsed(1))./(runAgesOfUsed(end)-runAgesOfUsed(1));

    %------ Calculate the sedrates for each pair of ages
    age_diffs = diff(runAgesOfUsed); %y
    dep_diffs = diff(depthOfUsed)'; %cm
    SRs = dep_diffs./age_diffs; %cm/y
    if S.normWithRunAve
        normSRs = SRs./aveSR_run;
    else
        normSRs = SRs./(scenario_meanSR(i_sce)./1000); %the ./1000 is a conversion from cm/kyr to cm/y to be consistent with units of SRs
    end

    %Add normSRs to vector to count them (with their weighting)
    weightingNormaliser = numruns; %Find normalising value based on number of runs
    nSRinfo = [normSRs; dep_diffs; age_diffs]; %Set up nSR info (nSR counts, depth differences, age differences)

    %Store all nSR info in a single array, with NaNs separating info from
    %different runs (instead of NaN in age-diff row, the first age used
    %is stored, so that the whole time history can be created from
    %this)
    if ix == 1
        nSRcounts = [[NaN; min(depthOfUsed); min(runAgesOfUsed)], nSRinfo];
        agediffs = age_diffs;
    else
        nSRcounts = cat(2, nSRcounts,[[NaN; min(depthOfUsed); min(runAgesOfUsed)], nSRinfo]);
        agediffs = cat(2, agediffs,age_diffs);
    end
end

%% Now need to create some illustrative plots

NaNlogi = isnan(nSRcounts(1,:));
NaNind = find(NaNlogi) ;
depths = cumsum(nSRcounts(2, NaNind(1):NaNind(2)-1));
ages = cumsum(nSRcounts(3, NaNind(1):NaNind(2)-1))./1000;

figure;
hold on;
plot(depths, ages, 'LineStyle', 'none', 'Marker', '+')
