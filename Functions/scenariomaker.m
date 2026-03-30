function [scenarios, chosen_labelsSplit, minus_labels] = scenariomaker(dup_depth_LabIDs, genrev_labIDs, label)
% scenariomaker  Enumerate all age-depth scenarios from a set of
%                mutually-exclusive date pairs.
%
% Given one or more pairs (or larger groups) of lab IDs where exactly one
% member must be selected per group, this function generates every possible
% combination and returns each combination as a "scenario" — a string
% vector of the lab IDs that would be included in that age-depth model run.
%
% Two sources of mutually-exclusive groups are supported and can be passed
% simultaneously:
%   - Doubly-dated depths (dup_depth_LabIDs): depths with more than one
%     date, where only one date can appear in any single scenario.
%   - Generic reversal pairs (genrev_labIDs): date pairs identified by
%     scenariosDealWithReversals as causing an age reversal, where one of
%     the pair must be removed.
%
% INPUTS
%   dup_depth_LabIDs - (cell array, 1×d) Each cell holds a string vector
%                      of lab IDs at one doubly-dated depth. Pass [] if
%                      there are no doubly-dated depths.
%   genrev_labIDs    - (cell array, 1×r) Each cell holds a 2×1 string
%                      vector of the two lab IDs in one reversal pair.
%                      Pass [] if called from scenariosDDD only.
%   label            - (string vector) Full list of lab IDs for this core
%                      (or scenario), used to identify the dates that are
%                      not involved in any DDD or reversal pair.
%
% OUTPUTS
%   scenarios          - (cell array, s×1) Each cell is a string vector of
%                        the lab IDs included in that scenario.
%   chosen_labelsSplit - (cell array, s×1) For each scenario, the lab IDs
%                        that were "chosen" (one per DDD/reversal group).
%   minus_labels       - (string vector) The lab IDs that are not part of
%                        any DDD or reversal group; these appear in every
%                        scenario unchanged.
%
% See also: scenariosDDD, scenariosDealWithReversals

%% Combine both sources of mutually-exclusive groups
labIDs_cell = [dup_depth_LabIDs genrev_labIDs];

%% Find the labels that are not in any exclusive group
% These appear unchanged in every scenario and can be computed once.
ind          = contains(string(label), vertcat(labIDs_cell{1,:}));
minus_labels = label(~ind);

%% Enumerate every combination (one choice per group)
combs      = combinations(labIDs_cell{:});
combs_size = size(combs);
scenarios          = cell(combs_size(1), 1);
chosen_labelsSplit = cell(combs_size(1), 1);

for x = 1:combs_size(1)
    chosen_labelsT      = combs(x,:);
    chosen_labelsJoined = join(chosen_labelsT{:,:}, '_');    % flatten table row to string
    chosen_labelsSplit{x} = split(chosen_labelsJoined, '_'); % split back to string vector
    scenarios{x}          = [minus_labels; chosen_labelsSplit{x}];
end

end
