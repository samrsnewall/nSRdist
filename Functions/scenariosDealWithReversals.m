function [SS, newScenIndicator, IDpairs, agediffV] = scenariosDealWithReversals(SS, depth_cm, ageprobAll, calAge, label, corename, duplicated_depths, IDpairs, agediffV, S, plotfigs)
% scenariosDealWithReversals  Check each scenario for age reversals and
%                             split any that contain them into sub-scenarios.
%
% Iterates through the current set of scenarios in SS. For each scenario
% not yet confirmed reversal-free (SS.CFR == 0), it calls scenariopdfNorm
% to compute the inverse SR PDF and detect reversals. The function processes
% scenarios one at a time:
%
%   - If a scenario has no reversals it is marked as confirmed (CFR = 1)
%     and its computed outputs are stored in SS.
%   - If a scenario has reversals, new sub-scenarios are generated (via
%     scenariomaker) and inserted into SS in place of the original. The
%     function then returns immediately with newScenIndicator = 1, so that
%     the caller (oneCoreScenarios) loops and re-checks from scratch.
%   - If more than 12 reversals are found the core is considered
%     unanalysable: SS is cleared and newScenIndicator = 0 is returned.
%
% Scenarios already confirmed as reversal-free (SS.CFR == 1) are passed
% through unchanged.
%
% INPUTS
%   SS               - scenarioStruct (see oneCoreScenarios). Fields:
%                        .scenarios     Cell array of LabID string vectors
%                        .CFR           Logical vector (1 = reversal-free)
%                        .chosenLabels  Cell array of chosen label vectors
%                        .invSRvals     Cell array of inverse-SR value vectors
%                        .invSRprobs    Cell array of inverse-SR probability vectors
%                        .aveSR        Numeric vector of ave SRs (cm/kyr)
%                        .numdatepairs  Numeric vector of date-pair counts
%                        .ageModes      Cell array of age-mode vectors
%                        .lengthSed     Numeric vector of sediment lengths (cm)
%                        .MSI_byage     Numeric vector of age-based MSI values
%                        .MSI_bydepth   Numeric vector of depth-based MSI values
%   depth_cm         - (numeric vector) Depth of each date (cm)
%   ageprobAll       - (matrix) Calibrated age probabilities for all dates
%   calAge           - (numeric vector) Calibrated age axis
%   label            - (string vector) Lab IDs after filtering
%   corename         - (string) Core identifier, used in warning messages
%   duplicated_depths - (numeric vector) Depths with more than one date
%   IDpairs          - (string array) Accumulator for date-pair IDs across
%                      calls (passed through to scenariopdfNorm)
%   agediffV         - (cell) Accumulator for age differences across calls
%   S                - (struct) Settings struct
%   plotfigs         - (logical) Whether to produce diagnostic plots
%
% OUTPUTS
%   SS               - Updated scenarioStruct
%   newScenIndicator - 1 if scenarios were modified (caller must loop again);
%                      0 if all scenarios are now confirmed reversal-free
%                      or if the core was abandoned due to too many reversals
%   IDpairs          - Updated date-pair ID accumulator
%   agediffV         - Updated age-difference accumulator
%
% See also: oneCoreScenarios, scenariopdfNorm, scenariomaker, removeDuplicateScenarios

%% Unpack scenarioStruct fields into local variables
scenarios           = SS.scenarios;
scenariosCFR        = SS.CFR;
chosenLabels        = SS.chosenLabels;
scenario_invSRvals  = SS.invSRvals;
scenario_invSRprobs = SS.invSRprobs;
scenario_aveSR     = SS.aveSR;
numdatepairs        = SS.numdatepairs;
ageModes            = SS.ageModes;
lengthSed           = SS.lengthSed;
MSI_byage           = SS.MSI_byage;
MSI_bydepth         = SS.MSI_bydepth;

%% Initialise output arrays (sized to match the current scenario count)
numScenarios         = length(scenarios);
scenario_invSRvals2  = cell(numScenarios, 1);
scenario_invSRprobs2 = cell(numScenarios, 1);
scenario_labels2     = cell(numScenarios, 1);
scenario_aveSR2     = nan(numScenarios, 1);
numdatepairs2        = nan(numScenarios, 1);
ageModes2            = cell(numScenarios, 1);
lengthSed2           = nan(numScenarios, 1);
lengthAge2           = nan(numScenarios, 1); %#ok<NASGU>
MSI_bydepth2         = nan(numScenarios, 1);
MSI_byage2           = nan(numScenarios, 1);
scenariosNew         = cell(numScenarios, 1);
scenariosNewCFR      = zeros(numScenarios, 1);
chosenLabels2        = cell(numScenarios, 1);

%% Process each scenario
for j = 1:numScenarios
    if scenariosCFR(j) == 0
        date_bool            = ismember(string(label), scenarios{j});
        date_is              = find(date_bool == 1);
        scenario_labels2{j}  = label(date_bool == 1);

        %Choose ages in this scenario
        ageprob = ageprobAll(:,date_is);

        %Check to make sure there aren't duplicated labels in the scenario. If
        %there are remove one
        [uniqueLabels, uniqueIdx] = unique(scenario_labels2{j});
        if length(uniqueIdx) ~= length(scenario_labels2{j})
            scenario_labels2{j} = uniqueLabels;
        end

        [scenario_invSRvals2{j}, scenario_invSRprobs2{j}, scenario_aveSR2(j),...
            reversalpairs, numdatepairs2(j), ageModes2{j}, lengthSed2(j),...
            lengthAge2(j), MSI_byage2(j), MSI_bydepth2(j), IDpairs, agediffV] =...
            scenariopdfNorm(depth_cm, date_is, label, ...
            ageprob, calAge, IDpairs, agediffV, S, plotfigs);

        if sum(size(ageModes2{j})) == 0
            a = 1; %#ok<NASGU>
        end

        %Find number of reversals
        numreversals = sum(reversalpairs);

        %------------------------------------------------------------------
        % Case 1: too many reversals — abandon this core entirely
        %------------------------------------------------------------------
        if numreversals > 12
            disp("Too many scenarios in core" + corename + ": " + num2str(numreversals))
            SS.scenarios    = cell([]);
            SS.CFR          = [];
            SS.chosenLabels = cell([]);
            SS.invSRvals    = cell([]);
            SS.invSRprobs   = cell([]);
            SS.aveSR       = nan(0,1);
            SS.numdatepairs = nan(0,1);
            SS.ageModes     = [];
            SS.lengthSed    = nan(0,1);
            SS.MSI_byage    = nan(0,1);
            SS.MSI_bydepth  = nan(0,1);
            newScenIndicator = 0;
            return
        end

        %------------------------------------------------------------------
        % Case 2: reversals exist — build new sub-scenarios and return
        %------------------------------------------------------------------
        if numreversals ~= 0
            % Note down the labIDs of the dates for which there is a reversal
            rev_labIDs = strings(numreversals, 2);
            for iii = 1:numreversals
                iv  = find(reversalpairs == 1);
                rev_labIDs(iii,1) = scenario_labels2{j}(iv(iii));
                rev_labIDs(iii,2) = scenario_labels2{j}(iv(iii)+1);
            end

            %Find out if any of the reversal pairings are related to
            %duplicated depths
            uniqueDepthLabels = label(~ismember(depth_cm, duplicated_depths));
            if ~isempty(duplicated_depths)
                indrevmin     = ismember(string(rev_labIDs), uniqueDepthLabels);
                OneDateDDDrev = sum(indrevmin, 2) == 1;

                if ismember(1, OneDateDDDrev)
                    %One date in the reversal pair is from a doubly-dated
                    %depth: remove the other date from the scenario
                    rev_labIDs     = rev_labIDs(OneDateDDDrev, :);
                    date2removeLog = ~ismember(rev_labIDs, chosenLabels{j});
                    date2remove    = rev_labIDs(date2removeLog);
                    dates2keepLog  = ~ismember(scenarios{j}, date2remove);
                    scenarios{j}   = scenarios{j}(dates2keepLog);

                    SS.scenarios    = scenarios;
                    SS.CFR          = [scenariosNewCFR(1:j-1); scenariosCFR(j:end)];
                    SS.chosenLabels = chosenLabels;
                    SS.invSRvals    = scenario_invSRvals2;
                    SS.invSRprobs   = scenario_invSRprobs2;
                    SS.aveSR       = scenario_aveSR2;
                    SS.numdatepairs = numdatepairs2;
                    SS.ageModes     = ageModes2;
                    SS.lengthSed    = lengthSed2;
                    SS.MSI_byage    = MSI_byage2;
                    SS.MSI_bydepth  = MSI_bydepth2;
                    newScenIndicator = 1;
                    return
                end
                %If any reversals have both dates not from a duplicated
                %depth (or both from duplicated depths), fall through to
                %the generic handler below.
            end

            %Generic case: split each reversal pair with choose-one /
            %leave-one logic to generate new sub-scenarios
            genrev_labIDs = cell(1, numreversals);
            for vi = 1:sum(reversalpairs)
                genrev_labIDs{1,vi} = rev_labIDs(vi,:)';
            end
            [newscenariosNR, ~, ~]     = scenariomaker([], genrev_labIDs, scenario_labels2{j});
            numNew                     = length(newscenariosNR);
            chosenLabelsNR(1:numNew,1) = chosenLabels(j);

            SS.scenarios    = [scenariosNew(1:j-1);    scenarios(j+1:end);    newscenariosNR];
            SS.CFR          = [scenariosNewCFR(1:j-1); scenariosCFR(j+1:end); zeros(size(newscenariosNR))];
            SS.chosenLabels = [chosenLabels2(1:j-1);   chosenLabels(j+1:end); chosenLabelsNR];
            SS.invSRvals    = [scenario_invSRvals(1:j-1);  scenario_invSRvals2(j+1:end);  cell(numNew,1)];
            SS.invSRprobs   = [scenario_invSRprobs(1:j-1); scenario_invSRprobs2(j+1:end); cell(numNew,1)];
            SS.aveSR       = [scenario_aveSR(1:j-1);  scenario_aveSR2(j+1:end);  nan(numNew,1)];
            SS.numdatepairs = [numdatepairs(1:j-1);     numdatepairs2(j+1:end);      nan(numNew,1)];
            SS.ageModes     = [ageModes(1:j-1);         ageModes2(j+1:end);          cell(numNew,1)];
            SS.lengthSed    = [lengthSed(1:j-1);        lengthSed2(j+1:end);         nan(numNew,1)];
            SS.MSI_byage    = [MSI_byage(1:j-1);        MSI_byage2(j+1:end);         nan(numNew,1)];
            SS.MSI_bydepth  = [MSI_bydepth(1:j-1);     MSI_bydepth2(j+1:end);       nan(numNew,1)];
            newScenIndicator = 1;
            return
        end

        %------------------------------------------------------------------
        % Case 3 (no reversals): mark scenario as confirmed and store results
        %------------------------------------------------------------------
        scenariosNew(j)    = scenarios(j);
        scenariosNewCFR(j) = 1;
        chosenLabels2(j)   = chosenLabels(j);

    else
        %Already confirmed reversal-free — copy all information straight through
        scenariosNew(j)      = scenarios(j);
        scenariosNewCFR(j)   = 1;
        ageModes2(j)         = ageModes(j);
        chosenLabels2(j)     = chosenLabels(j);
        scenario_invSRvals2{j}  = scenario_invSRvals{j};
        scenario_invSRprobs2{j} = scenario_invSRprobs{j};
        scenario_aveSR2(j)  = scenario_aveSR(j);
        numdatepairs2(j)     = numdatepairs(j);
        lengthSed2(j)        = lengthSed(j);
        MSI_bydepth2(j)      = MSI_bydepth(j);
        MSI_byage2(j)        = MSI_byage(j);
    end

    newScenIndicator = 0;
end

%% All scenarios confirmed — pack results back into output struct
SS.scenarios    = scenariosNew;
SS.CFR          = scenariosNewCFR;
SS.chosenLabels = chosenLabels2;
SS.invSRvals    = scenario_invSRvals2;
SS.invSRprobs   = scenario_invSRprobs2;
SS.aveSR       = scenario_aveSR2;
SS.numdatepairs = numdatepairs2;
SS.ageModes     = ageModes2;
SS.lengthSed    = lengthSed2;
SS.MSI_byage    = MSI_byage2;
SS.MSI_bydepth  = MSI_bydepth2;

end
