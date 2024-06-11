function [core_transnums, core_CSE2x, nSRcounts, agediffs] = oneCoreTMRestrict(corename, scenarios, LabIDs, incDepths, excLabIDs, excDepths, minAgeDiff, calcTM)
%% Return empties if core not used
%Return Nans if core not used
if isempty(scenarios) 
    core_transnums = nan(3,3);
    core_CSE2x = nan(3,1);
    nSRcounts = [];
    agediffs = [];
    return
end

%% Read in Radiocarbon Data
%Read in some radiocarbon data from a net cdf file
[age, depth_cm, error, label] = getDataWA(corename);

%% Filter ages
[age, depth_cm, error, label, emptybreak1, emptybreak2] = filtering(age, depth_cm, error, label, LabIDs, incDepths, excLabIDs, excDepths);

if emptybreak1 == 1 || emptybreak2 == 1
    core_transnums = nan(3,3);
    core_CSE2x = nan(3,1);
    nSRcounts = [];
    agediffs = [];
    return
end
    
%% Create fake labels if needed
% If there are no LabIDs, create temporary LabIDs so that scenarios can be made
if sum(contains(label, "NaN"))~=0
    disp("Creating Fake Labels in OneCoreTM for " + corename)
    for i = find(contains(label, "NaN"))
        label(i) = "FakeLabel" + num2str(i);
    end
end

%% Check which scenarios to analyse
%Set up empty vectors
validScenariosBool = ones(length(scenarios),1);

% Check which scenarios are invalid (due to reversals etc).
for i_sce = 1:length(scenarios)
    if scenarios{i_sce} == "invalid"
        validScenariosBool(i_sce) = 0;
    end
end
%% Loop through scenarios using random sampling SR calculations

%Initiate empty variables
scenarioTM = zeros(3,3,length(scenarios));
sce_transnums = zeros(3,3,length(scenarios));
sce_CSE2x = zeros(3,1,length(scenarios));


% Loop through scenarios to calculate TM
for i_sce = 1:length(scenarios)

    %Create empties if scenario is invalid
    if scenarios{i_sce} == "invalid"
        scenarioTM(:,:,i_sce) = NaN;
        continue
    end
    
    %Find the dates in the given scenario
    date_bool = contains(string(label), scenarios{i_sce});
    date_is = find(date_bool == 1);

    %% Perform Calibrations
    ageprob = multiMatcal(age, error, date_is);

    %% Set up years vector
    m20_years = 0:55000;

    %% Run a number of random samples, sampling each age with positivity rule
    %Set number of random samples
    numruns = 1500;
    transnums_allruns = zeros(3,3,numruns);
    CSE2x_allruns = zeros(3,1,numruns);
    datesUsedStore = zeros(numruns, length(date_is));
    for runN = 1:numruns
        %------ Create a potential run of ages
        run_age = NaN(1,length(date_is));
        for i = 1:length(date_is)
            %Set up age probabilities so that every age is older than the
            %age of sample shallower than it
            if i == 1
                poss_ages = m20_years;
                age_probs = ageprob(:,i);
                run_age(i) = randsample(poss_ages, 1, true, age_probs);
            else 
                min_age = max(run_age(1:i-1));
                possibleAgesLog = m20_years>min_age;
                poss_ages = m20_years(possibleAgesLog);
                age_probs = ageprob(:,i);
                age_probs = age_probs(possibleAgesLog);
                age_probs = age_probs./sum(age_probs);
                run_age(i) = randsample(poss_ages, 1, true, age_probs);
                %%% RESTRICTION CODE
                %If next age is not more than x years greater than
                %previous accepted age, reject it and move to sample from next
                %radiocarbon date (x stored in input minAgeDiff)
                if run_age(i) - max(run_age(1:i-1)) < minAgeDiff
                    run_age(i) = NaN;
                end
            end
        end

        %% Find out which ages were used and which weren't
        datesUsed = ~isnan(run_age);
        depthOfUsed = depth_cm(datesUsed);
        runAgesOfUsed = run_age(datesUsed);
        datesUsedStore(i_sce,:) = datesUsed;

        %% Calculate mean SR for that potential run
        meanSR_run = (depthOfUsed(end)-depthOfUsed(1))./(runAgesOfUsed(end)-runAgesOfUsed(1));

        %% Calculate the sedrates for each pair of ages
        age_diffs = diff(runAgesOfUsed);
        dep_diffs = diff(depthOfUsed)';
        SRs = dep_diffs./age_diffs;
        normSRs = SRs./meanSR_run;
        %Add normSRs to vector to count them (with their weighting)
        weightingNormaliser = (sum(validScenariosBool).*numruns); %Find normalising value based on number of scenarios and number of runs
        nSRinfo = [normSRs; dep_diffs./weightingNormaliser; age_diffs./weightingNormaliser]; %Set up nSR info (nSR counts, depth differences, age differences)
        
        %Store all nSR info in a single array, with NaNs separating info from
        %different runs
        if runN == 1 && i_sce == 1
            nSRcounts = [nSRinfo, NaN(3,1)];
            agediffs = age_diffs;
        else
            nSRcounts = cat(2, nSRcounts,[nSRinfo, NaN(3,1)]);
            agediffs = cat(2, agediffs,age_diffs);
        end

        %% Calculate Transition Matrix if desired
        if calcTM == false
            continue
        end
        %---- Categorise normalised sed rates
        %Categories are Steady, Expansion, Contraction (S, E, C)
        char_categ = '';
        for i = 1:(length(runAgesOfUsed)-1)
            if normSRs(i) >= 0.922 && normSRs(i) <1.085
                char_categ(i) = 'S';

            elseif normSRs(i)>=1.085 && normSRs(i) < inf
                    char_categ(i) = 'E';

            elseif normSRs(i) >= 0 && normSRs(i) < 0.922
                    char_categ(i) = 'C';
            else %If the normSR is not one of those values it'll be a NaN, which signifies the end of a run of ages.
                char_categ(i) = 'b';

            end
        end

        %----- Convert this to a transition matrix
        %Find the number of times each transition occurs
        transnums = NaN(3,3,1);
        for jj = 1:1
            transnums(1,1,jj) = length(strfind(char_categ, "CC")); %Probability of Contraction to Contraction
            transnums(1,2,jj) = length(strfind(char_categ, "CS")); %Probability of Contraction to Steady
            transnums(1,3,jj) = length(strfind(char_categ, "CE")); %Probability of Contraction to Expansion
            transnums(2,1,jj) = length(strfind(char_categ, "SC")); %Probability of Steady to Contraction
            transnums(2,2,jj) = length(strfind(char_categ, "SS")); %Probability of Steady to Steady
            transnums(2,3,jj) = length(strfind(char_categ, "SE")); %Probability of Steady to Expansion
            transnums(3,1,jj) = length(strfind(char_categ, "EC")); %Probability of Expansion to Contraction
            transnums(3,2,jj) = length(strfind(char_categ, "ES")); %Probability of Expansion to Steady
            transnums(3,3,jj) = length(strfind(char_categ, "EE")); %Probability of Expansion to Expansion
        end

        %Find the number of times transitions occur from each state
        numCSE2x = zeros(3,1);
        numCSE2x(1) = length(strfind(char_categ(1:end-1), "C"));
        numCSE2x(2) = length(strfind(char_categ(1:end-1), "S"));
        numCSE2x(3) = length(strfind(char_categ(1:end-1), "E"));
        
        %Store these numbers for all runs
        transnums_allruns(:,:,runN) = transnums;
        CSE2x_allruns(:,:,runN) = numCSE2x;
    end
    %Sum over all runs
    if calcTM == false
        continue
    end
    sce_transnums(:,:,i_sce) = sum(transnums_allruns, 3);
    sce_CSE2x(:,:,i_sce) = sum(CSE2x_allruns, 3);
    scenarioTM(:,:,i_sce) = sum(transnums_allruns, 3)./sum(CSE2x_allruns, 3);
end
if calcTM == false
    core_transnums = nan(3,3);
    core_CSE2x = nan(3,1);
    return
end
core_transnums = sum(sce_transnums, 3)./length(scenarios);
core_CSE2x = sum(sce_CSE2x, 3)./length(scenarios);
end