function[scenarios, duplicated_depths, chosenLabels] = scenariosDDD(depth_cm, label, corename)
%Scenarios Doubly Dated Depths - This function deals with doubly dated
%depths by creating scenarios that each only choose one of the dates at
%each of those depths

%Assumes that all data are ordered by ascending depth
depth_diffs             = diff(depth_cm);                                   %Find difference between neighbouring depth values (data should be ordered by ascending depths)
duplicated_firstinds    = depth_diffs <= 0.00001;                           %Find where neighbouring values don't differ (some leeway given in case data aren't processed perfectly)
duplicated_depths       = unique(depth_cm([duplicated_firstinds' false]));  %Find which depths are duplicated

if ~isempty(duplicated_depths) %Run loop if there are duplicately dated depths

    %Initialise cell
    dup_depth_LabIDs = cell(1, length(duplicated_depths));

    %Find labels for each depth
    for idupdepth = 1:length(duplicated_depths)
        duplicated_allinds          = ismember(depth_cm, duplicated_depths(idupdepth)); %Find the indices of all duplicated depths (duplicated_firstind ony gives indice of first instance of duplicated depth)
        dup_depth_LabIDs{idupdepth} = label(duplicated_allinds); %Find labels that correspond to duplicated depths.
    end
    
    %Display the core name, the depths that have been dated multiple times,
    %and the labIDs of the dates for each depth.
    % disp("Core " + string(corename) + " has " + num2str(length(duplicated_depths)) + " doubly-dated depth(s).")
    % for n_dupdepth = 1:length(duplicated_depths)
    %     disp("Depth " + join(num2str(duplicated_depths(n_dupdepth))) + " has dates with the following labIDs " + join(dup_depth_LabIDs{n_dupdepth}))
    % end

    [scenarios, chosenLabels, ~] = scenariomaker(dup_depth_LabIDs, [], label);
else
    %If there are no duplicately dated depths, only 1 scenario, which is
    %all dates
    scenarios{1} = label;
    chosenLabels{1} = [];
end