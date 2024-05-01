function[ChosenMSPF, EM, ManE, NotCCR, emptybreak1, emptybreak2] = filteringForPlotting(age, depth, error, label, LabIDs, incDepths, excLabIDs, excDepths)

emptybreak1 = 0;
emptybreak2 = 0;

%% Choose only MSPF dates
%Filter the data so it's only those with MSPF
if strcmp(string(LabIDs), "all") %if string says all, do nothing
elseif ~isempty(LabIDs) % if string is not empty, check for depths to include
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
    logi1 = contains(strip(string(label)), strip(string(LabIDsSplitGood)));

    %fill excluded material (EM) outputs
    EM.age = age(~logi1);
    EM.depth = depth(~logi1);
    EM.error = error(~logi1);
    EM.label = label(~logi1);
elseif ~isempty(incDepths) %if there are depths to include include those

    %Choose what depths to include
    incDepthsStr = string(split(incDepths, ', '));
    incDepthsN = double(incDepthsStr);
    logi1 = ismember(depth, incDepthsN.*100);

    %fill excluded material (EM) outputs
    EM.age = age(~logi1);
    EM.depth = depth(~logi1);
    EM.error = error(~logi1);
    EM.label = label(~logi1);


else %if there is no information to include, signal break out of function
    emptybreak1 = 1; %return signal to break out of parent function
    return           %break out of this function
end

%% Remove manually excluded (ManE) dates

%Remove the LabIDs of ages that are clearly erroneous (large age
%reversal - hand picked)
if ~isempty(excLabIDs)
    logi3 = contains(string(label), string(split(excLabIDs, ', ')));
    ManE.age = age(logi3);
    ManE.depth = depth(logi3);
    ManE.error = error(logi3);
    ManE.label = label(logi3);
elseif ~isempty(excDepths)
    excDepthsStr = string(split(excDepths, ', '));
    excDepthsN = double(excDepthsStr);
    logi3 = ismember(depth, excDepthsN.*100);
    ManE.age = age(logi3);
    ManE.depth = depth(logi3);
    ManE.error = error(logi3);
    ManE.label = label(logi3);
else
    %if there are no dates to manually exclude, move on
end

%Display warning if both excLabIDs and excDepths have information
if ~isempty(excLabIDs) & ~isempty(excDepths)
    warning("Both excLabIDs and excDepths have information. Currently, if excLabIDs has information, excDepths will be ignored.")
end
%% Filter ages to fit within calibration curve
%Filter the data so it only includes those where the age is greater than
%1000 and below 42,000 yr C14 BP. (to ensure they fit within the limits of
%the calibration curve)
logi4 = find(age-error <= 1 & age+error >=42);
NotCCR.depth = depth(~logi4);
NotCCR.age = age(~logi4);
NotCCR.error = error(~logi4);
NotCCR.label = label(~logi4);

%% Filter cores to have at least 4 dates remaining (arbitrary choice... chosen to exclude cores with v low data)
%Filter out cores that have less than 4 dates remaining
if length(age)<4
    emptybreak2 = 1;
end

ChosenInd = logi1 & ~logi3 & ~ logi4;

%fill outputs
ChosenMSPF.age = age(ChosenInd);
ChosenMSPF.depth = depth(ChosenInd);
ChosenMSPF.error = error(ChosenInd);
ChosenMSPF.label = label(ChosenInd);
end