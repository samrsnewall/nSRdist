function[modenSRinfo, mediannSRinfo, nSRcounts, LinNSRcounts, meanSR] = nSRBchron(corename,LabIDs, incDepths, excLabIDs, excDepths, dataLoc, S)
% nSRBchron  Compute normalised sedimentation rates (nSR) for a single core
%            using Bchron Bayesian age-depth modelling.
%
% Manages the full Bchron workflow for one sediment core:
%   1. Checks whether Bchron results already exist. Uses them if they do
%   (unless S.BchronReDo is true) - creates them if not.
%   2.1. Reads radiocarbon data from the appropriate source (World Atlas or
%      Lin2014 dataset).
%   2.2. Optionally filters dates using the same criteria as the RSR methods.
%   2.3. Writes a Bchron-formatted input file and calls Bchron via Rscript
%   3. Loads the Bchron output and computes nSR four ways:
%        BMode   - using the mode   of the Bchron age ensemble
%        BMedian - using the median of the Bchron age ensemble
%        BSamp   - using each iteration of the Bchron ensemble individually
%        LinNSR  - using exact method of Lin2014       
%
% INPUTS
%   corename   - (string) Sediment core identifier (e.g. "RC13-228")
%   LabIDs     - (string array) Laboratory IDs for special-case filtering
%                (passed to filtering.m)
%   incDepths  - (numeric vector) Depths (cm) forced to be included during
%                filtering (passed to filtering.m)
%   excLabIDs  - (string array) Laboratory IDs to exclude during filtering
%                (passed to filtering.m)
%   excDepths  - (numeric vector) Depths (cm) to exclude during filtering
%                (passed to filtering.m)
%   dataLoc    - (string) Where to get Radiocarbon data from: "WA" (World Atlas) or
%                "Lin2014"
%   S          - (struct) Settings struct. Relevant fields:
%                  .sandboxPath            Root path of the repository
%                  .BchronFolderName       Subfolder name under BchronFolders/
%                  .BchronCalCurve         Calibration curve (e.g. "Marine20")
%                  .BchronReDo             (logical) Force re-run of Bchron
%                  .BchronFilter           (logical) Apply date filtering
%                  .BchronOutlier          Prior outlier probability per date
%                  .BchronReversalCriteria Rejection threshold (fraction):
%                                          dates flagged as outliers in more
%                                          than this fraction of runs are
%                                          excluded (e.g. 0.5 = 50%)
%                  .BchronDepthSpacing     Depth resolution (cm) for age
%                                          predictions
%                  .DeltaRError            Marine reservoir age correction
%                  .RscriptPath            Full path to the Rscript executable
%                  .minNumberOfAges        Minimum accepted ages required per
%                                          core or per individual run
%                  .normWithRunAve         (logical) If true, normalise each
%                                          BSamp run by its own mean SR; if
%                                          false, normalise all runs by the
%                                          overall mean SR from the mode solution
%
% OUTPUTS
%   modenSRinfo   - (3 x N) nSR matrix computed using mode ages (BMode)
%   mediannSRinfo - (3 x N) nSR matrix computed using median ages (BMedian)
%   nSRcounts     - (3 x M) nSR matrix from individual Bchron MCMC iterations (BSamp);
%                   blocks from successive runs are concatenated horizontally
%   LinNSRcounts  - (3 x N) nSR matrix computed following approach of Lin
%                   et al. (2014) code
%   meanSR        - (scalar) Overall mean sedimentation rate (cm/kyr) from
%                   the depth/age range of the median-age solution
%
% See also: calcData, getDataWA, getDatatxt, filtering, oneCoreRSR

%--- Locate Bchron output folders ---
BchronFolder = S.BchronFolderName;

%Set up directory to where the Bchron Data is
BchronFolderPath = fullfile(S.sandboxPath, "BchronFolders", BchronFolder);

% If Bchron Folder doesn't yet exist, make it
if ~isfolder(BchronFolderPath)
    mkdir(BchronFolderPath)
end

% Set up path to this core
coreDir = fullfile(S.sandboxPath, "BchronFolders", BchronFolder, "Outputs",corename);

% --- See if Bchron needs to be run: if so, run it ----

%If the BchronData Output folder for this core doesn't exist, run Bchronology
if ~isfolder(coreDir) || S.BchronReDo 

    %If Bchron Data "/Outputs" folder doesn't exist, create it
    outputsFolder = fullfile(S.sandboxPath, "BchronFolders", BchronFolder, "Outputs");
    if ~isfolder(outputsFolder)
        mkdir(outputsFolder)
    end

    %If Bchron Data inputs folder doesn't exist, create it
    inputsFolder = fullfile(S.sandboxPath, "BchronFolders", BchronFolder, "Inputs");
    if ~isfolder(inputsFolder)
        mkdir(inputsFolder)
    end

    %Set up name of Bchron input file
    c14input = fullfile(S.sandboxPath, "BchronFolders", BchronFolder, "Inputs",corename + "_radiocarbon.txt");
    
    %If the input file doesn't exist, create it
    if ~isfile(c14input) || S.BchronReDo
        %Get the radiocarbon data from our sources (either World Atlas or
        %Lin2014 dataset)
        if dataLoc == "WA"
            [age, depth_cm, error, label] = getDataWA(corename, S);
        elseif dataLoc == "Lin2014"
            [age, depth_cm, error, label] = getDatatxt(corename, S);
        end

        % If desired to filter, similar to RSR methods, apply filter
        if S.BchronFilter
            [age, depth_cm, error, label, emptybreak1, emptybreak2,~,~,~] = filtering(age, depth_cm, error, label, LabIDs, incDepths, excLabIDs, excDepths, corename, S);
            if emptybreak1 == 1 || emptybreak2 == 1
                modenSRinfo = [];
                mediannSRinfo = [];
                nSRcounts = [];
                LinNSRcounts = [];
                meanSR = NaN;
                return
            end
        end

        %Write a new tab delimited file in the correct format for Bchron to
        %read, using the radiocarbon data we have
        numages = length(age);
        outlierLikelihood = S.BchronOutlier;
        T = table(label, age.*1000, error.*1000, depth_cm, zeros(numages,1),...
            ones(numages,1).*outlierLikelihood, ones(numages,1).*0.001, ...
            ones(numages,1),'VariableNames', ...
            ["id", "Age", "Error", "Depth", "Thickness", "Outlier1", "Outlier2", "Type"]);
        writetable(T, c14input, "Delimiter", 'tab')

    end

    %Bring in which calibration curve will be used
    calCurve = S.BchronCalCurve;

    %Run Bchronology by sending a line to the terminal (which allows us
    %to run Rscript within Matlab), with the name of the folder to save it
    %in, the name of the core, and the calibration curve to be used, the
    %reservoir age error, and the depth spacing to use
    disp("Running Bchronology for core " + corename)
    cmnd = S.RscriptPath + " " +...
        fullfile(S.sandboxPath, "Functions", "runBchron.R")...
        + " " + S.sandboxPath + " " + fullfile("BchronFolders", BchronFolder) + " "...
        + corename + " " + calCurve + " " + S.DeltaRError + " " + S.BchronDepthSpacing;
    [status, output] = system(cmnd);
    
    %Double check whether the command was run smoothly or not
    if status == 0
        disp("|============================================================| 100%")
    else
        disp("R output below:")
        disp(output)
        error("Problem in Bchron")
    end
end

% --- Use Bchron output files to get NSR matrix ---

%Load radiocarbon data
rData = readtable(fullfile(coreDir, "inputData.txt"));

%load individual run data
thetaData = readmatrix(fullfile(coreDir, "theta.csv"), "NumHeaderLines",1);
phiData = readmatrix(fullfile(coreDir, "phi.csv"), "NumHeaderLines", 1);

%Get the thetaData for the depths where radiocarbon ages were taken
depthR_r2BDS = round(rData.Depth./S.BchronDepthSpacing).*S.BchronDepthSpacing; %This is the depth of each radiocarbon date, rounded to the BchronDepthSpacing used for predictPosition
predictPositions = depthR_r2BDS(1):S.BchronDepthSpacing:depthR_r2BDS(end);    %This is all the depths Bchron estimated ages
predictPositionsRLOG = ismembertol(predictPositions,depthR_r2BDS);             %this is the logical to get ages estimated by Bchron at radiocarbon date depths
predictPositionsR = predictPositions(predictPositionsRLOG);                 %This is depths on predictPositions that are radiocarbon dated
thetaDataR = thetaData(:, predictPositionsRLOG);                            %This is ages at those depths

%Get mode and median from thetaData (which holds the 1000 ages used in the
%1000 Bchron iterations for each radiocarbon dated depth - it's different than the radiocarbon calibration).

%%% NOTE THAT THE MODE calculated using MATLAB is often different to the
%%% MODE using Bchron. This is possible when the dataset is multimodal
%%% (which the thetaData typically are) - MATLAB will always choose the
%%% lowest mode, whilst R will (depending on your function, because there
%%% is no base R function to get the mode) get whichever Mode occurs first
%%% in the dataset. The mode is a bad choice for finding a single choice in
%%% a numerical dataset because it depends on the resolution of the data
%%% (Bchron output is likelihood per year) and there may be multiple modes.
%%% 
modeMAT = mode(thetaDataR);
medianMAT = median(thetaDataR);
pctRejMAT = sum(phiData)./10;

%Test if there are any multi-modal Bchron estimates
%[~,~,test] = mode(thetaData);
% if any(cellfun(@length, test) > 1)
%     disp("Core " + corename + " has a multimodal Bchronology estimate")
% end

%Get relevant information for calculating nSR with mode of Ages
modeAge = modeMAT';
medianAge = medianMAT';
pctReject = pctRejMAT';

%Remove rejected ages if desired (Lin2014 did not do this)

%Remove the age if the age was rejected more than some percentage of the
%time. Note, even in a normal looking radiocarbon profile,
%good ages are rejected up to 5% of the time.
removeRejectedAges = true;
if removeRejectedAges
    rejectLog = pctReject > S.BchronReversalCriteria*100;  %Find which ages were rejected more than the criterium set
    acceptedDepths = depthR_r2BDS(~rejectLog);                  %Find out what ages were accepted (using their depth info to index them) 
    depthsUsed = unique(acceptedDepths);         %Find out what depths had estimated ages we want to keep (unique and accepted depths)   
    rejectAgeLog = ~ismembertol(depthsUsed, predictPositionsR); %Find out what estimated ages we want to reject (ages are estimated at unique depths, keep the unique depths we accepted)
    modeAgeUsed = modeAge(~rejectAgeLog);       %keep ages not rejected
    medianAgeUsed = medianAge(~rejectAgeLog);   %keep ages not rejected
else
    depthsUsed   = depthR_r2BDS;
    modeAgeUsed  = modeAge;
    medianAgeUsed = medianAge;
end

%Handle doubly-dated depths
%First option is to just find the unique depths, because each unique depth
%has a unique and singular modeAge value.
[depthsUsed, ia, ~]  = unique(depthsUsed);
modeAgeUsed = modeAgeUsed(ia);
medianAgeUsed = medianAgeUsed(ia);

%highlight problem if less than 4 ages are used
if length(depthsUsed) < S.minNumberOfAges
    warning(corename + " has less than " + num2str(S.minNumberOfAges) + " accepted ages after Bchron")
end

%Calculate nSR with mode ***************************
%calculate meanSR
meanSR = (depthsUsed(end)-depthsUsed(1))./(modeAgeUsed(end) - modeAgeUsed(1));
SRs = diff(depthsUsed)./diff(modeAgeUsed);
nSRs = SRs./meanSR;
depthdiffs = diff(depthsUsed);
agediffs = diff(modeAgeUsed);

%Store all nSR info in one matrix, the standard set up I use
modenSRinfo = [NaN, nSRs'; depthsUsed(1), depthdiffs'; modeAgeUsed(1), agediffs'];

%Calculate nSR with median ages **************************
meanSR = (depthsUsed(end)-depthsUsed(1))./(medianAgeUsed(end) - medianAgeUsed(1));
SRs = diff(depthsUsed)./diff(medianAgeUsed);
nSRs = SRs./meanSR;
agediffs = diff(medianAgeUsed);

%Store all nSR info in one matrix, the standard set up I use
mediannSRinfo = [NaN, nSRs'; depthsUsed(1), depthdiffs'; medianAgeUsed(1), agediffs'];

%%% Calculate nSR with probability of ages method *********************

% Find out which runs have less than 4 ages accepted
runs_with_less_than_x_ages = sum((size(phiData, 2)-sum(phiData,2))<S.minNumberOfAges);
if runs_with_less_than_x_ages > size(thetaDataR,1)/5
warning("Core " + corename + " has " + num2str(runs_with_less_than_x_ages) + "runs with less than " + num2str(S.minNumberOfAges) + " ages" )
end
%Go through each run, calculate the nSR info, and store it
for i = 1:size(thetaDataR,1)
    %Find out which ages were not rejected in this run (this also deals
    %with doubly dated depths, by finding unique depths and using thetaRs,
    %which only estimate one age per depth)
    keepAges = ~logical(phiData(i,:));

    %See if this run had enough ages (commented out because this reduces
    %the number of Bchron runs used to below 1000)
    % if sum(keepAges)<S.minNumberOfAges                
    %     continue
    % end

    %Get depths of radiocarbon ages not rejected in this run
    rundepthspDD = depthR_r2BDS(keepAges)'; 
    rundepths = unique(rundepthspDD);       

    %Get ages of those depths for this run
    keepThetaRs = ismembertol(rundepths, predictPositionsR);
    runages = thetaDataR(i,keepThetaRs);
    
    %Calculate SRs, depthdiffs, agediffs
    runSRs = diff(rundepths)./diff(runages);
    rundepthdiffs = diff(rundepths);
    runagediffs = diff(runages);

    %Calculate nSR
    if S.normWithRunAve
        runmeanSR = (rundepths(end)-rundepths(1))./(runages(end)-runages(1));
        runnSR = runSRs./runmeanSR;
    else
        runnSR = runSRs./meanSR;
    end

    %Store info in one block using the 3-row nSR matrix format
    nSRinfo = [NaN, runnSR; rundepths(1), rundepthdiffs; runages(1), runagediffs];
    if ~exist('nSRcounts','var')
        nSRcounts = nSRinfo;
    else
        nSRcounts = [nSRcounts, nSRinfo]; %#ok<AGROW>
    end
end

%%% Calculate NSR with lin2014 method %%%%%%%%%%%%%%%%%%%%%%%%%%%%

if S.Lin2014Method
% ------Following code reproduces data that Lin2014 code loaded in ------

% Compute median posterior age at each predictPosition
medAge = median(thetaData, 1);  % median ages (in yr) of each predictPosition

% Interpolate to find depth at each 1 kyr step of median age
age_targets = ceil(min(medAge)):floor(max(medAge));   % find the relevant 1kyr time steps
[medAge_unique, idx] = unique(medAge, 'stable'); %make sure medAge is unique - in theory it is, but because it is estimated by sampling, there can be some duplication, very certain this is not a serious problem
depth_at_age = interp1(medAge_unique, predictPositions(idx), age_targets, 'linear'); %find depths at those 1kyr time steps

% read in median calibrated ages
calT = readtable(fullfile(coreDir, "calAgesMedian.csv"));
medCalAge = calT.pct50; %median ages (yr) of each calibrated radiocarbon age (not influenced by Bchron)

% ------ Following code is exact replica of code used by Lin2014 ---------

% Bchron age models frequently don't include first and last 14c dates
% due to even 1-kyr increments of output
% Add these ages back in
if rData.Depth(1) < min(depth_at_age)
    depth_at_age=[rData.Depth(1), depth_at_age];
    age_targets =[medCalAge(1), age_targets];
end
if rData.Depth(end) > max(depth_at_age)
    depth_at_age=[depth_at_age, rData.Depth(end)];
    age_targets=[age_targets, medCalAge(end)];
end

%Bchron model ages for depth of each 14c measurement
Bage_14c=interp1(depth_at_age,age_targets,rData.Depth);

%Calculate sed rates (NaN indicates duplicate 14C at a particular depth)
cores_sr=diff(rData.Depth)./diff(Bage_14c);

% mean sed rate for core
cores_meansr=(max(rData.Depth)-min(rData.Depth))/(max(Bage_14c)-min(Bage_14c));

% Divide by mean sed rate to convert to sed rate ratio
sed_rate_ratio=cores_sr/cores_meansr;
ind1=find(~isnan(sed_rate_ratio));

% Duration of sed rate
depint=diff(rData.Depth);
dur=diff(Bage_14c); 
ind2=find(dur~=0); %SN Note: ind1 and ind2 will be identical
LinNSRcounts=[NaN, sed_rate_ratio(ind1)'; rData.Depth(1), depint(ind2)';  Bage_14c(1), dur(ind2)'];
else
    LinNSRcounts = [];
end