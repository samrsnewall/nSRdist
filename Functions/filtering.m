function[age, depth, error, label, emptybreak1, emptybreak2, EM, ManE, NotCCR, ageGapTooHigh] = filtering(age, depth, error, label, LabIDs, incDepths, excLabIDs, excDepths, corename, S)
% filtering  Filter radiocarbon dates for a single sediment core.
%
% Applies four successive filters to a core's radiocarbon dataset:
%   1. Material selection — keep only dates matching the specified LabIDs
%      (or depths in incDepths), typically used to retain MSPF dates only.
%   2. Manual exclusion — remove dates flagged by excLabIDs or excDepths
%      (e.g. known age reversals, non-planktonic material).
%   3. Calibration curve range — remove dates whose age ± error falls
%      outside the range S.c14AgeLim (ensures MatCal can calibrate them).
%   4. Minimum date count — sets emptybreak2 if fewer than S.minNumberOfAges
%      dates remain after the above filters.
%
% INPUTS
%   age       - (numeric vector) Radiocarbon ages (14C kyr BP)
%   depth     - (numeric vector) Sample depths (cm)
%   error     - (numeric vector) Radiocarbon age 1-sigma errors (14C kyr)
%   label     - (string vector) Laboratory IDs for each date
%   LabIDs    - (string) Comma-separated lab IDs to keep, or "all" to
%               retain all dates regardless of material type
%   incDepths - (string) Comma-separated depths (m) to keep; used as an
%               alternative to LabIDs when lab IDs are unavailable
%   excLabIDs - (string) Comma-separated lab IDs to exclude manually
%   excDepths - (string or numeric) Depths (m) to exclude manually
%   corename  - (string) Core identifier, used in warning messages
%   S         - (struct) Settings struct. Relevant fields:
%                 .c14AgeLim      [min, max] radiocarbon age limits (14C yr)
%                 .minNumberOfAges  Minimum number of dates required
%
% OUTPUTS
%   age, depth, error, label  - Filtered versions of the input arrays
%   emptybreak1  - (logical) 1 if no dates matched the material filter
%                  (LabIDs/incDepths empty or unmatched); caller should
%                  skip this core
%   emptybreak2  - (logical) 1 if fewer than S.minNumberOfAges dates
%                  remain after filtering; caller should skip this core
%   EM           - Struct of dates excluded by the material filter
%   ManE         - Struct of dates excluded by manual exclusion
%   NotCCR       - Struct of dates excluded for being outside the
%                  calibration curve range
%   ageGapTooHigh - Struct of dates involved in age gaps > 5 kyr
%                   (diagnostic only; not used to exclude dates)
%
% See also: calcData, oneCoreScenarios, oneCoreRSR, nSRBchron

emptybreak1 = 0; % This will signal if core has no data chosen (1 = no data)
emptybreak2 = 0; % This will signal if core has less than S.minNumberOfAges (1 = less than...)

%% Choose dates desired (i.e. planktonic foraminifera dates only)
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
        warning("Some of the LabIDs listed to be chosen in the COPYcorechoices_MSPF file for core" + corename + " do not match with the LabIDs read in from the WA2022")
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

logi31 = false(length(age), 1);
logi32 = false(length(age), 1);

%Remove the LabIDs of ages that are clearly erroneous (large age
%reversal - manually chosen)

%specified by LabID
if ~isempty(excLabIDs) & ~isnan(excLabIDs)
    logi31 = contains(string(label), string(split(excLabIDs, ', ')));
end

%specified by depth
if ~isempty(excDepths) & ~isnan(excDepths)
    if isa(excDepths, 'char')
        excDepthsStr = string(split(excDepths, ', '));
        excDepthsN = double(excDepthsStr);
    elseif isa(excDepths, 'double')
        excDepthsN = excDepths;
    else
        disp("excDepths data is not being used for a core because of data type - see filtering function")
    end
    logi32 = ismembertol(depth, excDepthsN.*100, 0.00001);
end

%Combine the two logicals and set these into Manually Excluded category
logi3 = logi31 ==1 |logi32 ==1;
ManE.age = age(logi3);
ManE.depth = depth(logi3);
ManE.error = error(logi3);
ManE.label = label(logi3);

%% Filter ages to fit within calibration curve
%Filter the data so it only includes those where the age is greater than
%1000 and below 42,000 yr C14 BP. (to ensure they fit within the limits of
%the calibration curve)
logi4 = age-error >= S.c14AgeLim(1) & age+error <=S.c14AgeLim(2);
NotCCR.depth = depth(~logi4);
NotCCR.age = age(~logi4);
NotCCR.error = error(~logi4);
NotCCR.label = label(~logi4);

%% Choose cores that fit all criteria
ChosenLogi = logi1 & ~logi3 & logi4;

%fill outputs
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

%% Break the core at place of large age gap and see if portion can be used?

%% Filter cores to have at least x dates remaining (arbitrary choice... chosen to exclude cores with v low data)
% %Filter out cores that have less than x dates remaining
if length(age)<S.minNumberOfAges
    emptybreak2 = 1;
end


end