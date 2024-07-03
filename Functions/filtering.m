function[age, depth, error, label, emptybreak1, emptybreak2, EM, ManE, NotCCR, ageGapTooHigh] = filtering(age, depth, error, label, LabIDs, incDepths, excLabIDs, excDepths)

emptybreak1 = 0;
emptybreak2 = 0;

%% Choose only MSPF dates
%Filter the data so it's only those with MSPF
if strcmp(string(LabIDs), "all") %if string says all, no dates are excluded for their material
    logi1 = true(length(age), 1);
    EM.age = age(~logi1);
    EM.depth = depth(~logi1);
    EM.error = error(~logi1);
    EM.label = label(~logi1);
elseif ~isempty(LabIDs) & ~isnan(LabIDs)  % if string is not empty, check for depths to include
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
    logi1 = ismember(strip(string(label)), strip(string(LabIDsSplitGood)));
    if sum(logi1)~=length(LabIDsSplitGood)
        warning("Some of the LabIDs listed to be chosen in the COPYcorechoices_MSPF file do not match with the LabIDs read in from the WA2022")
    end

    %fill excluded material (EM) outputs
    EM.age = age(~logi1);
    EM.depth = depth(~logi1);
    EM.error = error(~logi1);
    EM.label = label(~logi1);
elseif ~isempty(incDepths) & ~isnan(incDepths)  %if there are depths to include include those

    %Choose what depths to include
    incDepthsStr = string(split(incDepths, ', '));
    incDepthsN = double(incDepthsStr);
    logi1 = ismembertol(depth, incDepthsN.*100, 1e-6); %tolerance of 1e-6 due to low precision in core GEOFARK12

    %fill excluded material (EM) outputs
    EM.age = age(~logi1);
    EM.depth = depth(~logi1);
    EM.error = error(~logi1);
    EM.label = label(~logi1);

else %if there is no information to include, signal break out of function
    emptybreak1 = 1; %return signal to break out of parent function
    logi1 = false(size(age));
    EM.age = age;
    EM.depth = depth;
    EM.error = error;
    EM.label = label;
end

%% Remove manually excluded (ManE) dates

%Remove the LabIDs of ages that are clearly erroneous (large age
%reversal - hand picked)
if ~isempty(excLabIDs) & ~isnan(excLabIDs)
    logi3 = contains(string(label), string(split(excLabIDs, ', ')));
    ManE.age = age(logi3);
    ManE.depth = depth(logi3);
    ManE.error = error(logi3);
    ManE.label = label(logi3);
elseif ~isempty(excDepths) & ~isnan(excDepths)
    if isa(excDepths, 'char')
        excDepthsStr = string(split(excDepths, ', '));
        excDepthsN = double(excDepthsStr);
    elseif isa(excDepths, 'double')
        excDepthsN = excDepths;
    else
        disp("excDepths data is not being used for a core because of data type - see filtering function")
    end
    logi3 = ismember(depth, excDepthsN.*100);
    ManE.age = age(logi3);
    ManE.depth = depth(logi3);
    ManE.error = error(logi3);
    ManE.label = label(logi3);


else
    logi3 = false(length(age), 1);
    ManE.age = age(logi3);
    ManE.depth = depth(logi3);
    ManE.error = error(logi3);
    ManE.label = label(logi3);
    %if there are no dates to manually exclude, move on
end

%Display warning if both excLabIDs and excDepths have information
if ~isempty(excLabIDs) & ~isnan(excDepths)
    warning("Both excLabIDs and excDepths have information. Currently, if excLabIDs has information, excDepths will be ignored.")
end
%% Filter ages to fit within calibration curve
%Filter the data so it only includes those where the age is greater than
%1000 and below 42,000 yr C14 BP. (to ensure they fit within the limits of
%the calibration curve)
logi4 = age-error >= 1 & age+error <=42;
NotCCR.depth = depth(~logi4);
NotCCR.age = age(~logi4);
NotCCR.error = error(~logi4);
NotCCR.label = label(~logi4);



%% Choose cores that fit all criteria
ChosenLogi = logi1 & ~logi3 & logi4;

%fill outputs
% ChosenMSPF.age = age(ChosenLogi);
% ChosenMSPF.depth = depth(ChosenLogi);
% ChosenMSPF.error = error(ChosenLogi);
% ChosenMSPF.label = label(ChosenLogi);
age = age(ChosenLogi);
depth = depth(ChosenLogi);
error = error(ChosenLogi);
label = label(ChosenLogi);
%% Note where there are age differences greater than 5kyr (using mean of uncalibrated radiocarbon age)
agediffs = diff(age); %calculate differences between means of uncalibrated radiocarbon age
logi5 = agediffs > 5; %find locations where differences are greater than 5kyr
logi6 = false(size(age)); %initiate a logical (default = false) to denote which ages are involved in the differences greater than desired value
if sum(logi5) ~= 0
    for i = 1:length(agediffs)
        if logi5(i) == 1
            logi6(i) = 1;
            logi6(i+1) = 1;
        end
    end
end
ageGapTooHigh.depth = depth(logi6);
ageGapTooHigh.age = age(logi6);
ageGapTooHigh.error = error(logi6);
ageGapTooHigh.label = label(logi6);

%% Break the core at place of large age gap and see if portion can be used

%% Filter cores to have at least 4 dates remaining (arbitrary choice... chosen to exclude cores with v low data)
%Filter out cores that have less than 4 dates remaining
if length(age)<4
    emptybreak2 = 1;
end


end