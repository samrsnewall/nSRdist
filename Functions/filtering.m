function[age, depth_cm, error, label, emptybreak1, emptybreak2] = filtering(age, depth_cm, error, label, LabIDs, incDepths, excLabIDs, excDepths)

emptybreak1 = 0;
emptybreak2 = 0;

%% Choose only MSPF dates
%Filter the data so it's only those with MSPF
if strcmp(string(LabIDs), "all") %if string says all, do nothing
elseif isempty(LabIDs) % if string is empty, check for depths to include
    if isempty(incDepths) %if no depths and no LabIDs to include, break out of function
        emptybreak1 = 1; %return signal to break out of function
        return
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

    end
else %if string has LabIDs, find out which LabIDs they are by comparing with label
    
    %Get Good Lab IDs
    LabIDsSplit_Raw = split(LabIDs, ', ');
    
    %Remove problematic whitespace that isn't removed by strip (char 65279
    %causes many problems)
    test4ProblematicWhitespace = contains(LabIDsSplit_Raw, char(65279));
    if sum(test4ProblematicWhitespace) > 0
        disp("Removing problematic whitespace")
        LabIDsSplit_noProblemWhitespace = cell(length(LabIDsSplit_Raw),1);
        for i = 1:length(LabIDsSplit_Raw)
            LabIDcharProb = LabIDsSplit_Raw{i};
            problemIndex = ismember(LabIDcharProb, char(65279));
            LabIDcharFix = LabIDcharProb(~problemIndex);
            LabIDsSplit_noProblemWhitespace{i} = LabIDcharFix;
        end
        LabIDsSplitGood = LabIDsSplit_noProblemWhitespace;
    else
    LabIDsSplitGood = LabIDsSplit_Raw;
    end
    
    %Get indeces of labels that are in LabIDs field
    ind = contains(strip(string(label)), strip(string(LabIDsSplitGood)));
   
    %fill outputs with these indeces
    age = age(ind);
    depth_cm = depth_cm(ind);
    error = error(ind);
    label = label(ind);
end

%% Remove manually excluded dates

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

%% Filter ages to fit within calibration curve
%Filter the data so it only includes those where the age is greater than
%1000 and below 42,000 yr C14 BP. (to ensure they fit within the limits of
%the calibration curve)
ind = find(age-error >= 1 & age+error <=42);
depth_cm = depth_cm(ind);
age = age(ind);
error = error(ind);
label = label(ind);

%% Filter cores to have at least 4 dates remaining (arbitrary choice... chosen to exclude cores with v low data)
%Filter out cores that have less than 4 dates remaining
if length(age)<4
emptybreak2 = 1;
end

end