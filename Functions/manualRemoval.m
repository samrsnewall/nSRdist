function[age, depth_cm, error, label] = manualRemoval(age,depth_cm, error, label, excLabIDs, excDepths)
%Remove the LabIDs of ages that are clearly erroneous (large age
%reversal - hand picked)
if isempty(excLabIDs)
else
    ind = contains(string(label), string(split(excLabIDs, ', ')));
    ind = ind==0;
    age = age(ind);
    depth_cm = depth_cm(ind);
    error = error(ind);
    label = label(ind);
end

if isempty(excDepths)
else
    excDepthsStr = string(split(excDepths, ', '));
    excDepthsN = double(excDepthsStr);
    ind = ismember(depth_cm, excDepthsN.*100);
    ind = ind==0;
    age = age(ind);
    depth_cm = depth_cm(ind);
    error = error(ind);
    label = label(ind);
end
end
