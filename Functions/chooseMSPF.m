function[age, depth_cm, error, label, emptybreak] = chooseMSPF(age, depth_cm, error, label, LabIDs, incDepths)

%Filter the data so it's only those with MSPF
if strcmp(string(LabIDs), "all") %if string says all, do nothing
elseif isempty(LabIDs) % if string is empty, check for depths to include
    if isempty(incDepths) %if no depths and no LabIDs to include, break out of function
        emptybreak = 1; %return signal to break out of function
    else 
        %Choose what depths to include
        incDepthsStr = string(split(incDepths, ', '));
        incDepthsN = double(incDepthsStr);
        ind = ismember(depth_cm, incDepthsN.*100);

        %fill outputs
        age = age(ind);
        depth_cm = depth_cm(ind);
        error = error(ind);
        label = label(ind);
        emptybreak = 0;
    end
else %if string has LabIDs, find out which LabIDs they are by comparing with label
    ind = contains(string(label), string(split(LabIDs, ', ')));

    %fill outputs
    age = age(ind);
    depth_cm = depth_cm(ind);
    error = error(ind);
    label = label(ind);
    emptybreak = 0;
end