function[scenariosNew, scenariosNewCFR, chosenLabels2, scenario_invSRvals, scenario_invSRprobs, scenario_meanSR, reversalpairs, numDatePairs, ageModes, lengthSed,  newScenIndicator, MSI_byage, MSI_bydepth] = scenariosDealWithReversals(scenarios, scenariosCFR, chosenLabels, depth_cm, ageprobAll, calAge, label, corename, duplicated_depths, S, plotfigs)

%Initiate cells to hold results from each scenario
numScenarios        = length(scenarios);
scenario_invSRvals  = cell(numScenarios, 1);
scenario_invSRprobs = cell(numScenarios, 1);
scenario_labels     = cell(numScenarios, 1);
scenario_meanSR     = nan(numScenarios, 1);
numDatePairs        = nan(numScenarios, 1);
ageModes            = cell(numScenarios, 1);
lengthSed           = nan(numScenarios, 1);
lengthAge           = nan(numScenarios, 1);
MSI_bydepth         = nan(numScenarios, 1);
MSI_byage           = nan(numScenarios, 1);
scenariosNew        = cell(numScenarios, 1);
scenariosNewCFR     = zeros(numScenarios, 1);
chosenLabels2       = cell(numScenarios, 1);

%For each scenario, find the possible invSRvals and the probs
for j = 1:numScenarios
    if scenariosCFR(j) == 0
        date_bool           = ismember(string(label), scenarios{j});
        date_is             = find(date_bool == 1);
        scenario_labels{j}  = label(date_bool == 1);

        %Choose ages in this scenario
        ageprob = ageprobAll(:,date_is);

        %Check to make sure there aren't duplicated labels in the scenario. If
        %there are remove one
        [uniqueLabels, uniqueIdx] = unique(scenario_labels{j});
        if length(uniqueIdx) ~= length(scenario_labels{j})
            scenario_labels{j} = uniqueLabels;
        end

        [scenario_invSRvals{j}, scenario_invSRprobs{j}, scenario_meanSR(j),...
            reversalpairs, numDatePairs(j), ageModes{j}, lengthSed(j),...
            lengthAge(j), MSI_byage(j), MSI_bydepth(j)] =...
            scenariopdfNorm(depth_cm, date_is, ageprob, calAge, S, plotfigs);

        %Find number of reversals
        numreversals = sum(reversalpairs);

        if numreversals > 12
            %Will be too many scenarios, do not analyse core
            disp("Too many scenarios in core" + corename + ": " + num2str(numreversals))
            numScenarios = 0;
            scenario_invSRvals  = cell(numScenarios, 1);
            scenario_invSRprobs = cell(numScenarios, 1);
            % scenario_labels     = cell(numScenarios, 1);
            scenario_meanSR     = nan(numScenarios, 1);
            numDatePairs        = nan(numScenarios, 1);
            ageModes            = [];
            lengthSed           = nan(numScenarios, 1);
            % lengthAge           = nan(numScenarios, 1);
            MSI_bydepth         = nan(numScenarios, 1);
            MSI_byage           = nan(numScenarios, 1);
            scenariosNew        = cell([]);
            scenariosNewCFR     = [];
            chosenLabels2       = cell([]);
            newScenIndicator = 0;
            break
        end

        %If reversals exist, run through process to create new scenarios,
        %removing those reversals
        if numreversals ~=0
            % Note down the labIDs of the dates for which there is a reversal
            rev_labIDs = strings(numreversals,2);                          %Initiate
            for iii = 1:numreversals
                %find lab IDs for each reversal, each row represents a reversal
                %pairing
                iv  = find(reversalpairs == 1);
                rev_labIDs(iii,1) = scenario_labels{j}(iv(iii));
                rev_labIDs(iii,2) = scenario_labels{j}(iv(iii)+1);
            end

            %Find out if any of the reversal pairings are related to
            %duplicated depths
            uniqueDepthLabels        = label(~ismember(depth_cm, duplicated_depths));
            if ~isempty(duplicated_depths) %if there are duplicated depths
                indrevmin       = ismember(string(rev_labIDs), uniqueDepthLabels);   %find which labIDs are related to duplicated depths (0 = related, 1 = not related)
                OneDateDDDrev   = sum(indrevmin,2)==1;                          %Find which reversals have only one age from a duplicately dated depth

                if ismember(1,OneDateDDDrev) %If any reversals have one date from a duplicated depth then remove the date that isn't from the duplicated depth (for all such reversals)
                    rev_labIDs      = rev_labIDs(OneDateDDDrev, :);
                    date2removeLog  = ~ismember(rev_labIDs, chosenLabels{j});   %Find logical of date not from the duplicated depth, which will be removed
                    date2remove     = rev_labIDs(date2removeLog);               %Find label of date...
                    dates2keepLog   = ~ismember(scenarios{j}, date2remove);     %Find logicals of dates to keep from scenario
                    scenarios{j}    = scenarios{j}(dates2keepLog);              %Override the old scenario with the new scenario
                    scenariosNew    = scenarios;
                    scenariosNewCFR = scenariosCFR;
                    chosenLabels2   = chosenLabels;
                    newScenIndicator= 1;
                    break
                end
            %If any reversals have both dates not from a duplicated depth, update scenario with choose one leave one for those reversals
            %If any reversals have both dates from duplicated depths, update with choose one and leave one
            %Hence, if there are no reversals with only one date
            %part of a duplicated depth, then we want to update all
            %scenarios by leaving one choosing one for each
            %reversal, the same as how to treat if there are no
            %duplicate depths
            end

            %If there are no duplicate depths, all reversals must be generic, so update with choose one leave one
            genrev_labIDs = cell(1,numreversals);
            for vi = 1:sum(reversalpairs)
                genrev_labIDs{1,vi} = rev_labIDs(vi,:)';
            end
            [newscenariosNR, ~, ~]      = scenariomaker([], genrev_labIDs,scenario_labels{j});
            numNew                      = length(newscenariosNR);
            chosenLabelsNR(1:numNew,1)  = chosenLabels(j);
            scenariosNew                = [scenariosNew(1:j-1);    scenarios(j+1:end);    newscenariosNR];
            scenariosNewCFR             = [scenariosNewCFR(1:j-1); scenariosCFR(j+1:end); zeros(size(newscenariosNR))];
            chosenLabels2               = [chosenLabels2(1:j-1);   chosenLabels(j+1:end); chosenLabelsNR];
            newScenIndicator = 1;
            break
        end
    end

    %If the scenario shows no reversals, include it in the new scenario
    %list and label that it has been Checked For Reversals.
    scenariosNew(j)     = scenarios(j);
    scenariosNewCFR(j)  = 1;
    chosenLabels2(j)    = chosenLabels(j);
    newScenIndicator    = 0;
end