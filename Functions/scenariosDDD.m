function [scenarios, duplicated_depths, chosenLabels] = scenariosDDD(depth_cm, label, corename)
% scenariosDDD  Build age-depth scenarios for a core that has doubly-dated depths.
%
% Where two or more radiocarbon dates share exactly the same depth, only
% one can be used per age-depth model run. This function enumerates all
% possible "choose one per depth" combinations (via scenariomaker) so that
% each scenario contains exactly one date per depth. If there are no
% doubly-dated depths the function returns a single scenario containing all
% dates.
%
% A "doubly-dated depth" (DDD) is identified when two consecutive depth
% values differ by less than 0.00001 cm (a small tolerance to guard against
% floating-point imprecision in the input data). Only neighbouring duplicate
% pairs are caught by this tolerance check; data must be sorted by ascending
% depth beforehand.
%
% INPUTS
%   depth_cm  - (numeric vector, n×1) Depths of each date (cm), sorted
%               ascending
%   label     - (string vector, n×1) Lab IDs corresponding to each date
%   corename  - (string) Core identifier, used for diagnostic messages
%
% OUTPUTS
%   scenarios         - (cell array, s×1) Each cell holds a string vector of
%                       lab IDs that belong to that scenario
%   duplicated_depths - (numeric vector) Depths at which more than one date
%                       exists. Empty if there are no DDDs.
%   chosenLabels      - (cell array, s×1) For each scenario, the lab IDs
%                       that were "chosen" (i.e. one per DDD), so that the
%                       unchosen sibling can be identified later. Empty cell
%                       if there are no DDDs.
%
% See also: scenariomaker, oneCoreScenarios, scenariosDealWithReversals

%% Identify doubly-dated depths
% Assumes data are ordered by ascending depth.
depth_diffs          = diff(depth_cm);
duplicated_firstinds = depth_diffs <= 0.00001;
duplicated_depths    = unique(depth_cm([duplicated_firstinds' false]));

%% Build scenarios
if ~isempty(duplicated_depths)

    % Collect the lab IDs of every date at each duplicated depth.
    dup_depth_LabIDs = cell(1, length(duplicated_depths));
    for idupdepth = 1:length(duplicated_depths)
        duplicated_allinds         = ismember(depth_cm, duplicated_depths(idupdepth));
        dup_depth_LabIDs{idupdepth} = label(duplicated_allinds);
    end

    % Enumerate all combinations (one date per DDD) using scenariomaker.
    [scenarios, chosenLabels, ~] = scenariomaker(dup_depth_LabIDs, [], label);

else
    % No DDDs: only one scenario, consisting of all dates.
    scenarios{1}    = label;
    chosenLabels{1} = [];
end

end
