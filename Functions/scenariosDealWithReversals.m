function[scenariosNR, scenario_invSRvals, scenario_invSRprobs, scenario_meanSR, reversalpairs, numdatepairs, lengthsed,  newScenIndicator, res_byage, res_bydepth] = scenariosDealWithReversals(scenarios, depth_cm, age, error, label, corename, duplicated_depths, plotfigs)

%Initiate cells to hold results from each scenario
scenario_invSRvals = cell(1,length(scenarios));
scenario_invSRprobs = cell(1, length(scenarios));
scenario_labels = cell(1, length(scenarios));
scenario_meanSR = nan(1, length(scenarios));
numdatepairs = nan(1, length(scenarios));
lengthsed = nan(1, length(scenarios));
lengthage = nan(1, length(scenarios));
res_bydepth = nan(1, length(scenarios));
res_byage = nan(1, length(scenarios));
scenariosNR = cell([]);

% transprobs_scenarios = zeros(3,3,length(scenarios));
%For each scenario, find the possible invSRvals and the probs
for j = 1:length(scenarios)
    date_bool = contains(string(label), scenarios{j});
    date_is = find(date_bool == 1);
    scenario_labels{j} = label(date_bool == 1);
    [scenario_invSRvals{j}, scenario_invSRprobs{j}, scenario_meanSR(j), reversalpairs, numdatepairs(j), lengthsed(j), lengthage(j), res_byage(j), res_bydepth(j)] = scenariopdfNorm(depth_cm, age, error, scenario_labels{j}, date_is, plotfigs);
    %If a reversal arises, check whether it is related to any
    %duplicated depths or not
    numreversals = sum(reversalpairs);
    if numreversals ~=0
        disp("There are " + num2str(numreversals) + " reversals in scenario " + ...
            num2str(j) + " of core " + string(corename))
        rev_labIDs = strings(numreversals,2);
        for iii = 1:sum(reversalpairs) % note down the labIDs of the dates for which there is a reversal
            iv = find(reversalpairs == 1);
            %find lab IDs for each reversal, each row represents a reversal
            %pairing
            rev_labIDs(iii,1) = scenario_labels{j}(iv(iii));
            rev_labIDs(iii,2) = scenario_labels{j}(iv(iii)+1);
        end
        %Find out if any of the reversal pairings are related to
        %duplicated depths
        dup_depth_labels = label(ismember(depth_cm, duplicated_depths));
        minus_labels = label(~ismember(depth_cm, duplicated_depths));
        if ~isempty(duplicated_depths) %if there are duplicated depths
            indrevmin = contains(string(rev_labIDs), minus_labels); %find which labIDs are related to duplicated depths (0 = related, 1 = not related)
            generic_revs = sum(indrevmin,2)==2; %Find which reversals are both from ages not from duplicately dated depths
            if ismember(1, generic_revs) %If both dates from a pairing are not from a duplicated depth, break out of for loop and re-do scenarios, with choose-one-leave-one for the reversal
                genrev_labIDs = cell(1,sum(generic_revs));
                for vi = 1:sum(generic_revs)
                    gen_inds = find(generic_revs == 1);
                    genrev_labIDs{1,vi} = rev_labIDs(gen_inds(vi),:)';
                end
                cell_dup_depths = cell(1,1);
                cell_dup_depths{1} = dup_depth_labels;
                [scenariosNR, ~,~] = scenariomaker(cell_dup_depths, genrev_labIDs,label);
                newScenIndicator = 1;
                break
            else %If one of the dates is from a duplicated depth... Then what?
                disp("CODE NEEDS WORK IN FUNCTION scenariosDealWithReversals " + corename)

            end
        else %If there are no duplicate depths, all reversals must be generic
            genrev_labIDs = cell(1,sum(reversalpairs));
            for vi = 1:sum(reversalpairs)
                genrev_labIDs{1,vi} = rev_labIDs(vi,:)';
            end
            [twoscenariosNR, ~, ~] = scenariomaker([], genrev_labIDs,scenario_labels{j});
            scenariosNR = [scenariosNR; twoscenariosNR(1); twoscenariosNR(2)];
            newScenIndicator = 1;
            break
        end
    else
    scenariosNR = [scenariosNR; scenarios(j)];
    newScenIndicator = 0;
    end
end