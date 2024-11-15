function [nSRcounts, agediffs] = oneCoreTMRestrict(corename, scenarios, LabIDs, incDepths, excLabIDs, excDepths, scenario_meanSR, S, minAgeDiff)

%Check whether the core has previously been rejected
if ~isempty(scenarios) %If core has not been previously rejected...
    %------- Read in Radiocarbon Data
    %Read in some radiocarbon data from a net cdf file
    [age, depth_cm, error, label] = getDataWA(corename);

    %------- Filter ages
    [age, depth_cm, error, label, emptybreak1, emptybreak2] = filtering(age, depth_cm, error, label, LabIDs, incDepths, excLabIDs, excDepths, corename, S);

    if ~(emptybreak1 == 1 || emptybreak2 == 1) %If filtering of ages does not lead to core rejection
        %------ Create fake labels if needed
        % If there are no LabIDs, create temporary LabIDs so that scenarios can be made
        if sum(contains(label, "NaN"))~=0
            for i = find(contains(label, "NaN"))
                label(i) = "FakeLabel" + num2str(i);
            end
        end

        %Check to make sure there aren't duplicated labels in the scenario
        [uniqueLabels, uniqueIdx] = unique(label);
        if length(uniqueIdx) ~= length(label)
            label = uniqueLabels;
            for ilabel = 1:length(label)
                label(ilabel) = "xFakeLabel" + num2str(ilabel);
            end
        end

        %------- Check which scenarios to analyse
        %Set up empty vectors
        validScenariosBool = ones(length(scenarios),1);

        %Check which scenarios are invalid (due to reversals etc).
        for i_sce = 1:length(scenarios)
            if scenarios{i_sce} == "invalid"
                validScenariosBool(i_sce) = 0;
            end
        end

        %------- Loop through scenarios using random sampling SR calculations

        %Initiate empty variables
        scenarioTM = zeros(3,3,length(scenarios));
        % sce_transnums = zeros(3,3,length(scenarios));
        % sce_CSE2x = zeros(3,1,length(scenarios));

        % Loop through scenarios to calculate TM
        for i_sce = 1:length(scenarios)

            %Create empties if scenario is invalid
            if scenarios{i_sce} == "invalid"
                scenarioTM(:,:,i_sce) = NaN;
                continue
            end

            %Find the dates in the given scenario
            date_bool = ismember(string(label), scenarios{i_sce});
            date_is = find(date_bool == 1);
            

            % Perform Calibrations
            [ageprob, calAge] = multiMatcal(age, error, date_is, S);

            %get depths in this scenario
            depths_scenario = depth_cm(date_is);

            %------ Run a number of random samples, sampling each age with positivity rule
            %Set number of random samples
            numruns = 1000;
            %initialize
            % transnums_allruns = zeros(3,3,numruns);
            % CSE2x_allruns = zeros(3,1,numruns);
            datesUsedStore = zeros(numruns, length(date_is));
            for runN = 1:numruns
                %------ Create a potential run of ages
                run_age = NaN(1,length(date_is));
                for i = 1:length(date_is)
                    %Set up age probabilities so that every age is older than the
                    %age of sample shallower than it
                    if i == 1
                        poss_ages   = calAge;
                        age_probs   = ageprob(:,i);
                        run_age(i)  = randsample(poss_ages, 1, true, age_probs);
                    else
                        min_age = max(run_age(1:i-1));
                        possibleAgesLog = calAge>min_age;
                        poss_ages = calAge(possibleAgesLog);
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

                % Find out which ages were used and which weren't
                datesUsed = ~isnan(run_age);
                depthOfUsed = depths_scenario(datesUsed);
                runAgesOfUsed = run_age(datesUsed);
                datesUsedStore(runN,:) = datesUsed;

                % Calculate mean SR for that potential run
                meanSR_run = (depthOfUsed(end)-depthOfUsed(1))./(runAgesOfUsed(end)-runAgesOfUsed(1));

                %------ Calculate the sedrates for each pair of ages
                age_diffs = diff(runAgesOfUsed);
                dep_diffs = diff(depthOfUsed)';
                SRs = dep_diffs./age_diffs;
                if S.normWithRunMean
                normSRs = SRs./meanSR_run;
                else
                normSRs = SRs./(scenario_meanSR(i_sce)./1000);
                end
                %Add normSRs to vector to count them (with their weighting)
                weightingNormaliser = (sum(validScenariosBool).*numruns); %Find normalising value based on number of scenarios and number of runs
                nSRinfo = [normSRs; dep_diffs./weightingNormaliser; dep_diffs; age_diffs]; %Set up nSR info (nSR counts, weighting, depth differences, age differences)

                %Store all nSR info in a single array, with NaNs separating info from
                %different runs (instead of NaN in age-diff row, the first age used
                %is stored, so that the whole time history can be created from
                %this)
                if runN == 1 && i_sce == 1
                    nSRcounts = [[NaN; NaN; min(depthOfUsed); min(runAgesOfUsed)], nSRinfo];
                    agediffs = age_diffs;
                else
                    nSRcounts = cat(2, nSRcounts,[[NaN; NaN; min(depthOfUsed); min(runAgesOfUsed)], nSRinfo]);
                    agediffs = cat(2, agediffs,age_diffs);
                end

            end

            % [uniqueDatesUsed, ia, ic] = unique(datesUsedStore, 'rows');
            % [uDUoccurences, uDUindex] = groupcounts(ic);

        end

    else
        %Core was rejected after filtering through ages
        % core_transnums = nan(3,3);
        % core_CSE2x = nan(3,1);
        nSRcounts = [];
        agediffs = [];
        return
    end

else
    %Core was rejected previously
    % core_transnums = nan(3,3);
    % core_CSE2x = nan(3,1);
    nSRcounts = [];
    agediffs = [];
    return
end

end