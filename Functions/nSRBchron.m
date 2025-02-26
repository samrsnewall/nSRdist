function[modenSRinfo, mediannSRinfo, nSRcounts] = nSRBchron(corename, dataLoc, S)
%Find Bchron output data
BchronFolder = S.BchronFolderName;

%Set up directory to where the Bchron Data is
coreDir = fullfile(S.sandboxPath, BchronFolder, "Outputs",corename);

%Bring in which calibration curve will be used
calCurve = S.BchronCalCurve;

%If the BchronData Output folder doesn't exist, run Bchronology
if ~isfolder(coreDir) || S.BchronReDo 

    %If Bchron Data "/Outputs" folder doesn't exist, create it
    outputsFolder = fullfile(S.sandboxPath, BchronFolder, "Outputs");
    if ~isfolder(outputsFolder)
        mkdir(outputsFolder)
    end

    %If Bchron Data inputs folder doesn't exist, create it
    inputsFolder = fullfile(S.sandboxPath, BchronFolder, "Inputs");
    if ~isfolder(inputsFolder)
        mkdir(inputsFolder)
    end

    %Set up name of Bchron input file
    c14input = fullfile(S.sandboxPath, BchronFolder, "Inputs",corename + "_radiocarbon.txt");
    
    %If the input file doesn't exist, create it
    if ~isfile(c14input) || S.BchronReDo
        %Get the radiocarbon data from our sources (either World Atlas or
        %Lin2014 dataset)
        if dataLoc == "WA"
            [age, depth_cm, error, label] = getDataWA(corename, S);
        elseif dataLoc == "Lin2014"
            [age, depth_cm, error, label] = getDatatxt(corename, S);
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
        + " " + S.sandboxPath + " " + BchronFolder + " "...
        + corename + " " + calCurve + " " + S.DeltaRError;
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

%Get mode and median from thetaData (which holds the 1000 ages used in the
%1000 Bchron iterations for each radiocarbon date - it's different than the radiocarbon calibration).

%%% NOTE THAT THE MODE calculated using MATLAB is often different to the
%%% MODE using Bchron. This is possible when the dataset is multimodal
%%% (which the thetaData typically are) - MATLAB will always choose the
%%% lowest mode, whilst R will (depending on your function, because there
%%% is no base R function to get the mode) get whichever Mode occurs first
%%% in the dataset. The mode is a bad choice for finding a single choice in
%%% a numerical dataset because it depends on the resolution of the data
%%% (Bchron output is likelihood per year) and there may be multiple modes.
%%% 
modeMAT = mode(thetaData);
medianMAT = median(thetaData);
pctRejMAT = sum(phiData)./10;

%Test if there are any multi-modal Bchron estimates
[~,~,test] = mode(thetaData);
if any(cellfun(@length, test) > 1)
    disp("Core " + corename + " has a multimodal Bchronology estimate")
end

%Get relevant information for calculating nSR with mode of Ages
depths = rData.Depth;
modeAge = modeMAT';
medianAge = medianMAT';
pctReject = pctRejMAT';

%Remove rejected ages if desired (not sure how Lin2014 handled this)
%First option is to remove the mode of the age if the age was rejected more
%than 50% of the time. Note, even in a normal looking radiocarbon profile,
%good ages are rejected up to 5% of the time
removeRejectedAges = true;
if removeRejectedAges
    rejectLog = pctReject > 50;
    depthsUsed = depths(~rejectLog);
    modeAgeUsed = modeAge(~rejectLog);
    medianAgeUsed = medianAge(~rejectLog);
end

%Handle doubly-dated depths (not sure how Lin2014 handled this)
%First option is to just find the unique depths, because each unique depth
%has a unique and singular modeAge value.
[depthsUsed, ia, ~]  = unique(depthsUsed);
modeAgeUsed = modeAgeUsed(ia);
medianAgeUsed = medianAgeUsed(ia);

%Not sure how doubly dated depths will play through with the individual
%runs (as long as R is saving thetaPredict, not theta, then doubly dated depths will be given the same age)
% thetaData = thetaData(:,ia);
% phiData = phiData(:,ia);

%Calculate nSR with mode

%calculate meanSR
meanSR = (depthsUsed(end)-depthsUsed(1))./(modeAgeUsed(end) - modeAgeUsed(1));
SRs = diff(depthsUsed)./diff(modeAgeUsed);
nSRs = SRs./meanSR;
weights = diff(depthsUsed);
agediffs = diff(modeAgeUsed);

%Store all nSR info in one matrix, the standard set up I use
modenSRinfo = [NaN, nSRs'; NaN, weights'; depthsUsed(1), weights'; modeAgeUsed(1), agediffs'];

%Calculate nSR with median ages
meanSR = (depthsUsed(end)-depthsUsed(1))./(medianAgeUsed(end) - medianAgeUsed(1));
SRs = diff(depthsUsed)./diff(medianAgeUsed);
nSRs = SRs./meanSR;
weights = diff(depthsUsed);
agediffs = diff(medianAgeUsed);

%Store all nSR info in one matrix, the standard set up I use
mediannSRinfo = [NaN, nSRs'; NaN, weights'; depthsUsed(1), weights'; medianAgeUsed(1), agediffs'];

%%% Including probability of ages method
%Go through each run, calculate the nSR info, and store it
for i = 1:size(thetaData,1)
    %Find out which ages were not rejected in this run
    keepAges = ~logical(phiData(1,:));
    runagespDD = thetaData(i,keepAges);
    rundepthspDD = depths(keepAges)';

    %deal with doubly dated depths (will create duplicate runage and
    %rundepth values)
    runages = unique(runagespDD);
    rundepths = unique(rundepthspDD);
    
    %Calculate SRs, weights, agediffs
    runSRs = diff(rundepths)./diff(runages);
    runweights = diff(rundepths);
    runagediffs = diff(runages);

    %Calculate nSR
    if S.normWithRunMean
        runmeanSR = (rundepths(end)-rundepths(1))./(runages(end)-runages(1));
        runnSR = runSRs./runmeanSR;
    else
        runnSR = runSRs./meanSR;
    end

    %Store info in one vector, with my standardised format
    nSRinfo = [NaN, runnSR; NaN, runweights; rundepths(1), runweights; runages(1), runagediffs];
    if i == 1
        nSRcounts = nSRinfo;
    else
        nSRcounts = [nSRcounts, nSRinfo]; %#ok<AGROW>
    end
end