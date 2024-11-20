function[scenariosNew, scenariosNewCFR, chosenLabels2, scenario_invSRvals2, scenario_invSRprobs2, scenario_meanSR2, numdatepairs2, ageModes2, lengthSed2,  newScenIndicator, MSI_byage2, MSI_bydepth2, IDpairs, agediffV] = scenariosDealWithReversals(scenarios, scenariosCFR, chosenLabels, scenario_invSRvals, scenario_invSRprobs, scenario_meanSR, numdatepairs, ageModes, lengthSed, MSI_byage, MSI_bydepth, depth_cm, ageprobAll, calAge, label, corename, duplicated_depths, IDpairs, agediffV, S, plotfigs)

%Initiate cells to hold results from each scenario
numScenarios        = length(scenarios);
scenario_invSRvals2  = cell(numScenarios, 1);
scenario_invSRprobs2 = cell(numScenarios, 1);
scenario_labels2     = cell(numScenarios, 1);
scenario_meanSR2     = nan(numScenarios, 1);
numdatepairs2       = nan(numScenarios, 1);
ageModes2            = cell(numScenarios, 1);
lengthSed2           = nan(numScenarios, 1);
lengthAge2           = nan(numScenarios, 1);
MSI_bydepth2         = nan(numScenarios, 1);
MSI_byage2           = nan(numScenarios, 1);
scenariosNew        = cell(numScenarios, 1);
scenariosNewCFR     = zeros(numScenarios, 1);
chosenLabels2       = cell(numScenarios, 1);

%For each scenario, find the possible invSRvals and the probs
for j = 1:numScenarios
    if scenariosCFR(j) == 0
        date_bool           = ismember(string(label), scenarios{j});
        date_is             = find(date_bool == 1);
        scenario_labels2{j}  = label(date_bool == 1);

        %Choose ages in this scenario
        ageprob = ageprobAll(:,date_is);

        %Check to make sure there aren't duplicated labels in the scenario. If
        %there are remove one
        [uniqueLabels, uniqueIdx] = unique(scenario_labels2{j});
        if length(uniqueIdx) ~= length(scenario_labels2{j})
            scenario_labels2{j} = uniqueLabels;
        end

        [scenario_invSRvals2{j}, scenario_invSRprobs2{j}, scenario_meanSR2(j),...
            reversalpairs, numdatepairs2(j), ageModes{j}, lengthSed2(j),...
            lengthAge2(j), MSI_byage2(j), MSI_bydepth2(j), IDpairs, agediffV] =...
            scenariopdfNorm(depth_cm, date_is, label, ...
            ageprob, calAge, IDpairs, agediffV, S, plotfigs);

        if sum(size(ageModes{j})) == 0
            a = 1;
        end


        %Find number of reversals
        numreversals = sum(reversalpairs);

        if numreversals > 12
            %Will be too many scenarios, do not analyse core
            disp("Too many scenarios in core" + corename + ": " + num2str(numreversals))
            numScenarios = 0;
            scenario_invSRvals2  = cell(numScenarios, 1);
            scenario_invSRprobs2 = cell(numScenarios, 1);
            % scenario_labels     = cell(numScenarios, 1);
            scenario_meanSR2     = nan(numScenarios, 1);
            numdatepairs2        = nan(numScenarios, 1);
            ageModes2            = [];
            lengthSed2           = nan(numScenarios, 1);
            % lengthAge           = nan(numScenarios, 1);
            MSI_bydepth2         = nan(numScenarios, 1);
            MSI_byage2           = nan(numScenarios, 1);
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
                rev_labIDs(iii,1) = scenario_labels2{j}(iv(iii));
                rev_labIDs(iii,2) = scenario_labels2{j}(iv(iii)+1);
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
                    scenario_invSRvals2 = scenario_invSRvals;
                    scenario_invSRprobs2 = scenario_invSRvals;
                    scenario_meanSR2 = scenario_meanSR;
                    numdatepairs2 = numdatepairs;
                    ageModes2 = ageModes;
                    lengthSed2 = lengthSed;
                    MSI_bydepth2 = MSI_bydepth;
                    MSI_byage2 = MSI_bydepth;
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
            [newscenariosNR, ~, ~]      = scenariomaker([], genrev_labIDs,scenario_labels2{j});
            numNew                      = length(newscenariosNR);
            chosenLabelsNR(1:numNew,1)  = chosenLabels(j);
            scenariosNew                = [scenariosNew(1:j-1);    scenarios(j+1:end);    newscenariosNR];
            scenariosNewCFR             = [scenariosNewCFR(1:j-1); scenariosCFR(j+1:end); zeros(size(newscenariosNR))];
            chosenLabels2               = [chosenLabels2(1:j-1);   chosenLabels(j+1:end); chosenLabelsNR];
            scenario_invSRvals2         = [scenario_invSRvals(1:j-1); scenario_invSRvals2(j+1:end); cell(numNew,1)];
            scenario_invSRprobs2         = [scenario_invSRprobs(1:j-1); scenario_invSRprobs2(j+1:end); cell(numNew,1)];
            scenario_meanSR2         = [scenario_meanSR(1:j-1); scenario_meanSR2(j+1:end); nan(numNew,1)];
            numdatepairs2           = [numdatepairs(1:j-1); numdatepairs2(j+1:end); nan(numNew, 1)];
            ageModes2         = [ageModes(1:j-1); ageModes2(j+1:end); cell(numNew,1)];
            lengthSed2         = [lengthSed(1:j-1); lengthSed2(j+1:end); nan(numNew,1)];
            MSI_bydepth2         = [MSI_bydepth(1:j-1); MSI_bydepth2(j+1:end); nan(numNew,1)];
            MSI_byage2         = [MSI_byage(1:j-1); MSI_byage2(j+1:end); nan(numNew,1)];
            newScenIndicator = 1;
            break
        end
    

    %If the scenario shows no reversals, include it in the new scenario
    %list and label that it has been Checked For Reversals.
    scenariosNew(j)     = scenarios(j);
    scenariosNewCFR(j)  = 1;
    ageModes2(j)        = ageModes(j);
    chosenLabels2(j)    = chosenLabels(j);
    
    %Transfer over all other information
    else
    scenariosNew(j)     = scenarios(j);
    scenariosNewCFR(j)  = 1;
    ageModes2(j)        = ageModes(j);
    chosenLabels2(j)    = chosenLabels(j);
    scenario_invSRvals2{j} = scenario_invSRvals{j};
    scenario_invSRprobs2{j} = scenario_invSRprobs{j};
    scenario_meanSR2(j) = scenario_meanSR(j);
    numdatepairs2(j) = numdatepairs(j);
    lengthSed2(j) = lengthSed(j);
    MSI_bydepth2(j) = MSI_bydepth(j);
    MSI_byage2(j) = MSI_byage(j);
    end
    
    newScenIndicator    = 0;
end