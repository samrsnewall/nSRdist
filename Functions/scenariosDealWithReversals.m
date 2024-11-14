function[scenariosNew, scenariosNewCFR, chosenLabels2, scenario_invSRvals, scenario_invSRprobs, scenario_meanSR, reversalpairs, numDatePairs, ageModes, lengthSed,  newScenIndicator, MSI_byage, MSI_bydepth] = scenariosDealWithReversals(scenarios, scenariosCFR, chosenLabels, depth_cm, age, error, label, corename, duplicated_depths, S, plotfigs)

%Initiate cells to hold results from each scenario
numScenarios        = length(scenarios);
scenario_invSRvals  = cell(numScenarios, 1);
scenario_invSRprobs = cell(numScenarios, 1);
scenario_labels     = cell(numScenarios, 1);
scenario_meanSR     = nan(numScenarios, 1);
numDatePairs        = nan(numScenarios, 1);
ageModes            = [];
lengthSed           = nan(numScenarios, 1);
lengthAge           = nan(numScenarios, 1);
MSI_bydepth         = nan(numScenarios, 1);
MSI_byage           = nan(numScenarios, 1);
scenariosNew        = cell([]);
scenariosNewCFR     = [];
chosenLabels2       = cell([]);
% transprobs_scenarios = zeros(3,3,numscenarios);

%For each scenario, find the possible invSRvals and the probs
for j = 1:numScenarios
    if scenariosCFR(j) == 0;
        date_bool           = ismember(string(label), scenarios{j});
        date_is             = find(date_bool == 1);
        scenario_labels{j}  = label(date_bool == 1);

        %Check to make sure there aren't duplicated labels in the scenario. If
        %there are remove one
        [uniqueLabels, uniqueIdx] = unique(scenario_labels{j});
        if length(uniqueIdx) ~= length(scenario_labels{j});
            scenario_labels{j} = uniqueLabels;
        end

        [scenario_invSRvals{j}, scenario_invSRprobs{j}, scenario_meanSR(j),...
            reversalpairs, numDatePairs(j), ageModesScen, lengthSed(j),...
            lengthAge(j), MSI_byage(j), MSI_bydepth(j)] =...
            scenariopdfNorm(depth_cm, age, error, scenario_labels{j}, date_is,...
            S, plotfigs);

        %Find number of reversals
        numreversals = sum(reversalpairs);

        if numreversals > 12
            %Will be too many scenarios, do not analyse core
            disp("Too many scenarios in core" + corename + ": " + num2str(numreversals))
            numScenarios = 0;
            scenario_invSRvals  = cell(numScenarios, 1);
            scenario_invSRprobs = cell(numScenarios, 1);
            scenario_labels     = cell(numScenarios, 1);
            scenario_meanSR     = nan(numScenarios, 1);
            numDatePairs        = nan(numScenarios, 1);
            ageModes            = [];
            lengthSed           = nan(numScenarios, 1);
            lengthAge           = nan(numScenarios, 1);
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
            %Display how many reversals there are to command window.
            % disp("There are " + num2str(numreversals) + " reversals in scenario " + ...
            %     num2str(j) + " of core " + string(corename))

            rev_labIDs = strings(numreversals,2);

            for iii = 1:numreversals % note down the labIDs of the dates for which there is a reversal
                iv  = find(reversalpairs == 1);
                %find lab IDs for each reversal, each row represents a reversal
                %pairing
                rev_labIDs(iii,1) = scenario_labels{j}(iv(iii));
                rev_labIDs(iii,2) = scenario_labels{j}(iv(iii)+1);
            end
            %Find out if any of the reversal pairings are related to
            %duplicated depths
            dup_depth_labels    = label( ismember(depth_cm, duplicated_depths));
            uniqueDepthLabels        = label(~ismember(depth_cm, duplicated_depths));
            if ~isempty(duplicated_depths) %if there are duplicated depths
                indrevmin       = ismember(string(rev_labIDs), uniqueDepthLabels);   %find which labIDs are related to duplicated depths (0 = related, 1 = not related)
                generic_revs    = sum(indrevmin,2)==2;                          %Find which reversals are both from ages not from duplicately dated depths
                if ~ismember(0, generic_revs) %If both dates from a pairing are not from a duplicated depth, update scenario with choose one leave one, break out of loop
                    genrev_labIDs = cell(1,sum(generic_revs));
                    for vi = 1:sum(generic_revs)
                        gen_inds            = find(generic_revs == 1);
                        genrev_labIDs{1,vi} = rev_labIDs(gen_inds(vi),:)';
                    end
                    cell_dup_depths     = cell(1,1);
                    cell_dup_depths{1}  = dup_depth_labels;
                    [newscenariosNR, ~,~] = scenariomaker([], genrev_labIDs,scenario_labels{j});
                    chosenLabelsNR(1:length(newscenariosNR),1) = chosenLabels(j);
                    scenariosNew = [scenariosNew; scenarios(j+1:end); newscenariosNR];
                    scenariosNewCFR = [scenariosNewCFR; scenariosCFR(j+1:end); zeros(size(newscenariosNR))];
                    chosenLabels2 = [chosenLabels2; chosenLabels(1:1-j); chosenLabels(j+1:end); chosenLabelsNR];
                    newScenIndicator    = 1;

                    if length(scenariosNew) ~=length(scenariosNewCFR)
                        1+1
                    end
                    break
                else %If one of the dates is from a duplicated depth then remove the date that isn't from the duplicated depth
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
            else %If there are no duplicate depths, all reversals must be generic
                genrev_labIDs = cell(1,numreversals);
                for vi = 1:sum(reversalpairs)
                    genrev_labIDs{1,vi} = rev_labIDs(vi,:)';
                end
                [newscenariosNR, ~, ~] = scenariomaker([], genrev_labIDs,scenario_labels{j});
                chosenLabelsNR(1:length(newscenariosNR),1) = chosenLabels(j);
                    scenariosNew = [scenariosNew; scenarios(j+1:end); newscenariosNR];
                    scenariosNewCFR = [scenariosNewCFR; scenariosCFR(j+1:end); zeros(size(newscenariosNR))];
                    chosenLabels2 = [chosenLabels2; chosenLabels(1:1-j); chosenLabels(j+1:end); chosenLabelsNR];
                newScenIndicator = 1;

                break
            end
        else
            % scenariosNew = [scenariosNew; scenarios(j)];
            % scenariosNewCFR = [scenariosNewCFR; 1];
            % chosenLabels2 = [chosenLabels2; chosenLabels(j)];
            % ageModes = [ageModes;ageModesScen(~ismember(ageModesScen, ageModes))];
            % newScenIndicator = 1;
            % break
        end
    end
    scenariosNew = [scenariosNew; scenarios(j)];
    scenariosNewCFR = [scenariosNewCFR; 1];
    chosenLabels2 = [chosenLabels2; chosenLabels(j)];
    %ageModes = [ageModes;ageModesScen(~ismember(ageModesScen, ageModes))];
    newScenIndicator = 0;
end