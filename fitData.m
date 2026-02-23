%% Analyse a given set of results
%This script imports some results file (containing BMode, BSamp, etc...)
%results, and then fits distributions to them and performs statistical
%analysis using the chi-squared goodness of fit test.

%Add important paths
addpath('Functions')
addpath('Results')

%% Load the results
filepath =    "Results/dataT_All1_RLGtrue_DS0p05_Dec9";
load(filepath)
d.dataT = dataT;
d.S = S;
d.label = "AllCores, R200M20";
split_filepath = split(filepath, ".mat");
filename2save =  split_filepath(1)+ "_fit13Feb26_depthweight_400R.mat";

%% Set up fitting settings structure
fitS.Lin2014AgeFilter = [500 4000];
fitS.merge_small_dt = false;
fitS.weighting = "depth";    %Can be "depth", "age", or "none"
fitS.minMeanSR = 8;
fitS.minDepth = 1000;
fitS.maxLat = 40;

fitS.fitDists = true;
fitS.mlnReps = 3;           %For weighting with depth, best is 3
fitS.run_chi2gof = false;
fitS.chi2binN = 5;
fitS.dispChi2 = false;
fitS.chi2MinCountNum = 5;   %Minimum allowed counts in a bin (chisquared)
fitS.chi2MaxCountNum = 50;  %Maximum allowed counts in a bin (chisquared)
fitS.enforceBinSizeLimits = true; %(chisquared)
fitS.invXbinEdges = 0:0.1:15;
fitS.DeterministicRun.weightDP = 3;
fitS.DeterministicRun.weightInflator = 1;
fitS.resampleData = false;
fitS.OneRun.numruns = d.S.numruns;
if fitS.weighting == "depth"
    fitS.OneRun.weightRepDP = 3;       %For weighting with depth, best is 3
    fitS.OneRun.weightRepInflator = 1; %For weighting with depth, best is 4
    fitS.OneRun.MLNReps = 3;           %For weighting with depth, best is 3
elseif fitS.weighting == "age"
    fitS.OneRun.weightRepDP = 3;
    fitS.OneRun.weightRepInflator = 0.001;
    fitS.OneRun.MLNReps = 3;
else
    fitS.OneRun.weightRepDP = 1;       
    fitS.OneRun.weightRepInflator = 1; 
    fitS.OneRun.MLNReps = 3;           
end
fitS.divideLikelihood = true;
d.S1.fitS = fitS;


fitS.useParallelComp = false;
%% Fit Mix Log Norm to BIGMACS data for comparison reasons
%Load BIGMACS files
BMlognorm = readtable("BIGMACSdata/lognormal.txt");                  % (use x values currently used in BIGMACS)
BM.nSR.x = BMlognorm.Var1';
BM.nSR.px = BMlognorm.Var2;
BM.hist = readmatrix("BIGMACSdata/Lin2014_sedrateratio_cm_wo_NaN.txt");
BM.hist = BM.hist(:,4);
BM.TM = readmatrix("BIGMACSdata/transition_parameter.txt");
[MLN_BIGMACS, ~, BM.gmfit] = fitMixLogNorm(BM.hist, BM.nSR.x, 2, 3, 343);
[BM.lnSR.x, BM.lnSR.px] = px_to_pfx(BM.nSR.x, BM.nSR.px, @log);

%% Analyse new data
%Choose subset of cores
d.S1.chooseLog = d.dataT.meanSR > fitS.minMeanSR & d.dataT.depths > fitS.minDepth & abs(d.dataT.lats) < fitS.maxLat;
%Find number of cores
d.S1.numCores = sum(d.S1.chooseLog);

%% Fit MLN and inverse Gamma distribution to my Bchron Median data
%Set up specific fitting settings
fitS.Lin2014AgeFiltering = true;

%Create x vector that spans very wide range of nSR
d.logx = -6:0.01:6;
d.x = exp(d.logx);

%[d.S1.BMedian]= SRun_MLNandInvGam(d.dataT.bchronMedian, d.S1.chooseLog, 1, d.x, fitS);
[d.S1.BMedian]=ARfitdists(d.dataT.bchronMedian, d.x, d.S1.chooseLog, fitS.DeterministicRun.weightDP,fitS.DeterministicRun.weightInflator,1, fitS);
%% Fit dists to all BSamp samples
countDivisor = 1000;
[d.S1.BChAR] = ARfitdists(d.dataT.bchronProb, d.x, d.S1.chooseLog, 3, 1e-1, 1000, fitS);

%% Fit dists to individual Bchron Samplings
% %Fit many individual runs
fitS.Lin2014AgeFiltering = true;
numruns = fitS.OneRun.numruns;
fitS.dispChi2 = false;
rng(2)
[d.S1.BChIR] = SRun_MLNandInvGam(d.dataT.bchronProb, d.S1.chooseLog, numruns, d.x, fitS);

%% Fit dists to collected Newall Samplings
%weight replicator = 1e-1 keeps the mean and variance the same as NewIR,
%but is faster than weightReplicator = 1;
fitS.Lin2014AgeFiltering = false;
fitS.dispChi2 = false;
countDivisor = 400;
[d.S1.New0AR]    = ARfitdists(d.dataT.nSRcounts, d.x, d.S1.chooseLog, 3,1e-1, countDivisor, fitS); 
[d.S1.New500AR]  = ARfitdists(d.dataT.nSRcounts500, d.x, d.S1.chooseLog, 3, 1e-1,  countDivisor, fitS);
[d.S1.New1000AR] = ARfitdists(d.dataT.nSRcounts1000, d.x, d.S1.chooseLog, 3, 1e-1, countDivisor, fitS);

%% Fit dists to individual Newall Samplings
fitS.Lin2014AgeFiltering = false;
fitS.dispChi2 = false;
rng(2)
[d.S1.New0IR]    = SRun_MLNandInvGam(d.dataT.nSRcounts, d.S1.chooseLog, numruns, d.x, fitS);
[d.S1.New500IR]  = SRun_MLNandInvGam(d.dataT.nSRcounts500, d.S1.chooseLog, numruns, d.x, fitS);
[d.S1.New1000IR] = SRun_MLNandInvGam(d.dataT.nSRcounts1000, d.S1.chooseLog, numruns, d.x, fitS);

%% Calculate histogram counts in log Space for each individual run
bw           = 0.1; %bin width
logBinEdges  = -5:bw:5;          %lognSR
d.S1.BChIR.lnSRHistCounts = NaN(numruns, length(logBinEdges)-1);
d.S1.New0IR.lnSRHistCounts = NaN(numruns, length(logBinEdges)-1);
d.S1.New500IR.lnSRHistCounts = NaN(numruns, length(logBinEdges)-1);
d.S1.New1000IR.lnSRHistCounts = NaN(numruns, length(logBinEdges)-1);

for i = 1:numruns
    d.S1.BChIR.lnSRHistCounts(i,:)    = histcounts(log(d.S1.BChIR.weightedC{i}), 'BinEdges', logBinEdges)./d.S1.fitS.OneRun.weightRepInflator;
    d.S1.New0IR.lnSRHistCounts(i,:)   = histcounts(log(d.S1.New0IR.weightedC{i}), 'BinEdges', logBinEdges)./d.S1.fitS.OneRun.weightRepInflator;
    d.S1.New500IR.lnSRHistCounts(i,:) = histcounts(log(d.S1.New500IR.weightedC{i}), 'BinEdges', logBinEdges)./d.S1.fitS.OneRun.weightRepInflator;
    d.S1.New1000IR.lnSRHistCounts(i,:)= histcounts(log(d.S1.New1000IR.weightedC{i}), 'BinEdges', logBinEdges)./d.S1.fitS.OneRun.weightRepInflator;
end

%% Save data as file
if isfile(filename2save)
    disp("Warning, not saving file because it already exists. Please rename file")
else
    save(filename2save, "d", "-v7.3")
end