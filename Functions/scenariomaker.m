function[scenarios, chosen_labelsSplit, minus_labels] = scenariomaker(dup_depth_LabIDs, genrev_labIDs, label)

%Take cells arrays of pairs of LabIDs for which only one can be chosen for a given scenario
%Finds all possible scenarios with this data
labIDs_cell = [dup_depth_LabIDs genrev_labIDs];
combs = combinations(labIDs_cell{:}); %Create table of all possible pairings of dates to include
combs_size = size(combs);
scenarios = cell(combs_size(1),1); %Create holding space for each scenario (will have combinations of labIDs)
chosen_labelsSplit = cell(combs_size(1),1);
for x = 1:combs_size(1)
    all_labels = label; %bring in all lab ids for this core
    ind = contains(string(label), vertcat(labIDs_cell{1,:})); %find indices of labels that contain the lab IDs of duplicately dated depths
    minus_labels = all_labels(~ind); %Create a vector containing labels that are not from duplicately dated depths
    chosen_labelsT = combs(x,:); %Bring in the xth combination of labIDs (comes in a table)
    chosen_labelsJoined = join(chosen_labelsT{:,:}, '_'); %this step and next step removes it from a table
    chosen_labelsSplit{x} = split(chosen_labelsJoined, '_'); % see comment above
    chosen_labelsALL = [minus_labels; chosen_labelsSplit{x}]; % add labIDs from combination to those that aren't in duplicately dated depths
    scenarios{x} = chosen_labelsALL; % These labIDs create one scenario
end
end