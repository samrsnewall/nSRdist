function[modenSRinfo, nSRcounts] = nSRBchron(corename, S)

%Find Bchron output data, load mode Data
coreDir = fullfile("/Users/samnewall/Documents/MATLAB/nSRdist_code/Lin2014Cores/BchronOutputM09",corename);
modeData = readtable(fullfile(coreDir, "modeChron.csv"), "ReadVariableNames",true);
%load individual run data
thetaData = readmatrix(fullfile(coreDir, "theta.csv"), "NumHeaderLines",1);
phiData = readmatrix(fullfile(coreDir, "phi.csv"), "NumHeaderLines", 1);

%Get relevant information for calculating nSR with mode of Ages
depths = modeData.Var2;
modeAge = modeData.modeChron;
pctReject = modeData.rejectedAgePct;

%Handle doubly-dated depths (not sure how Lin2014 handled this)
[depths, ia, ic]  = unique(depths);
modeAge = modeAge(ia);
pctReject = pctReject(ia);
thetaData = thetaData(:,ia);
phiData = phiData(:,ia);

%Remove rejected ages if desired (not sure how Lin2014 handled this)
removeRejectedAges = true;
if removeRejectedAges
    rejectLog = pctReject > 50;
    depthsUsed = depths(~rejectLog);
    modeAgeUsed = modeAge(~rejectLog);
end

%calculate meanSR
meanSR = (depthsUsed(end)-depthsUsed(1))./(modeAgeUsed(end) - modeAgeUsed(1));
SRs = diff(depthsUsed)./diff(modeAgeUsed);
nSRs = SRs./meanSR;
weights = diff(depthsUsed);
agediffs = diff(modeAgeUsed);

modenSRinfo = [NaN, nSRs'; NaN, weights'; depthsUsed(1), weights'; modeAgeUsed(1), agediffs'];

%Go through each run, calculate the nSR info, and store it
for i = 1:size(thetaData,1)
    keepAges = ~logical(phiData(1,:));
    runages = thetaData(i,keepAges);
    rundepths = depths(keepAges)';

    runSRs = diff(rundepths)./diff(runages);
    runweights = diff(rundepths);
    runagediffs = diff(runages);
    if S.normWithRunMean
        runmeanSR = (rundepths(end)-rundeptphs(1))./(runages(end)-runages(1));
        runnSR = runSRs./runmeanSR;
    else
        runnSR = runSRs./meanSR;
    end

    nSRinfo = [NaN, runnSR; NaN, runweights; rundepths(1), runweights; runages(1), runagediffs];
    if i == 1
        nSRcounts = nSRinfo;
    else
        nSRcounts = [nSRcounts, nSRinfo];
    end
end