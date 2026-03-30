function[modenSRinfo, mediannSRinfo, nSRcounts, meanSR] = nSRBchron(corename,LabIDs, incDepths, excLabIDs, excDepths, dataLoc, S)
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
%   3. Loads the Bchron output and computes nSR three ways:
%        BMode   - using the mode   of the Bchron age ensemble
%        BMedian - using the median of the Bchron age ensemble
%        BSamp   - using each iteration of the Bchron ensemble individually
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
%   modenSRinfo   - (4 x N) nSR matrix computed using mode ages (BMode)
%   mediannSRinfo - (4 x N) nSR matrix computed using median ages (BMedian)
%   nSRcounts     - (4 x M) nSR matrix from individual Bchron runs (BSamp);
%                   blocks from successive runs are concatenated horizontally
%   meanSR        - (scalar) Overall mean sedimentation rate (cm/kyr) from
%                   the depth/age range of the median-age solution
%
% See also: calcData, getDataWA, getDatatxt, filtering, oneCoreRSR

%--- Locate Bchron output data
BchronFolder = S.BchronFolderName;

%Set up directory to where the Bchron Data is
coreDir = fullfile(S.sandboxPath, "BchronFolders", BchronFolder, "Outputs",corename);

%Bring in which calibration curve will be used
calCurve = S.BchronCalCurve;

%If the BchronData Output folder doesn't exist, run Bchronology
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

    %Run Bchronology by sending a line to the terminal (which will allow us
    %to run Rscript within Matlab), with the name of the folder to save it
    %in, the name of the core, and the calibration curve to be used
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

%Remove rejected ages if desired (not sure how Lin2014 handled this)
%First option is to remove the mode of the age if the age was rejected more
%than 50% of the time. Note, even in a normal looking radiocarbon profile,
%good ages are rejected up to 5% of the time
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

%Handle doubly-dated depths (not sure how Lin2014 handled this)
%First option is to just find the unique depths, because each unique depth
%has a unique and singular modeAge value.
[depthsUsed, ia, ~]  = unique(depthsUsed);
modeAgeUsed = modeAgeUsed(ia);
medianAgeUsed = medianAgeUsed(ia);

%highlight problem if less than 4 ages are used in bchron median
if length(depthsUsed) < S.minNumberOfAges
    warning(corename + " has less than " + num2str(S.minNumberOfAges) + " accepted ages after Bchron")
end

%For the individual runs, as long as R is saving thetaPredict, not theta, then doubly dated depths will be given the same age)

%Calculate nSR with mode ***************************
%calculate meanSR
meanSR = (depthsUsed(end)-depthsUsed(1))./(modeAgeUsed(end) - modeAgeUsed(1));
SRs = diff(depthsUsed)./diff(modeAgeUsed);
nSRs = SRs./meanSR;
depthdiffs = diff(depthsUsed);
agediffs = diff(modeAgeUsed);
weights = diff(depthsUsed);

%Store all nSR info in one matrix, the standard set up I use
modenSRinfo = [NaN, nSRs'; NaN, weights'; depthsUsed(1), depthdiffs'; modeAgeUsed(1), agediffs'];

%Calculate nSR with median ages **************************
meanSR = (depthsUsed(end)-depthsUsed(1))./(medianAgeUsed(end) - medianAgeUsed(1));
SRs = diff(depthsUsed)./diff(medianAgeUsed);
nSRs = SRs./meanSR;
weights = diff(depthsUsed);
agediffs = diff(medianAgeUsed);

%Store all nSR info in one matrix, the standard set up I use
mediannSRinfo = [NaN, nSRs'; NaN, weights'; depthsUsed(1), weights'; medianAgeUsed(1), agediffs'];

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
    if sum(keepAges)<S.minNumberOfAges
        continue
    end
    rundepthspDD = depthR_r2BDS(keepAges)';
    rundepths = unique(rundepthspDD);
    keepThetaRs = ismembertol(rundepths, predictPositionsR);
    runages = thetaDataR(i,keepThetaRs);
    
    %Calculate SRs, weights, agediffs
    runSRs = diff(rundepths)./diff(runages);
    rundepthdiffs = diff(rundepths);
    runweights = diff(rundepths);
    runagediffs = diff(runages);

    %Calculate nSR
    if S.normWithRunAve
        runmeanSR = (rundepths(end)-rundepths(1))./(runages(end)-runages(1));
        runnSR = runSRs./runmeanSR;
    else
        runnSR = runSRs./meanSR;
    end

    %Store info in one block using the standard 4-row nSR matrix format
    %(see function header for full format description)
    nSRinfo = [NaN, runnSR; NaN, runweights; rundepths(1), rundepthdiffs; runages(1), runagediffs];
    if ~exist('nSRcounts','var')
        nSRcounts = nSRinfo;
    else
        nSRcounts = [nSRcounts, nSRinfo]; %#ok<AGROW>
    end
end