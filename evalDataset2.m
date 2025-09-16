%% Analyse a given set of results
%This script imports some results file (containing BMode, BSamp, etc...)
%results, and then fits distributions to them and performs statistical
%analysis using the chi-squared goodness of fit test.

%Add important paths
addpath('../Functions')
addpath('../Results')

%% Load the results
filepath = "../Results/dataT_All1_RLGtrue_DS0p05_Jun2.mat";
load(filepath)
d.dataT = dataT;
d.S = S;
d.label = "AllCores, R200M20";
split_filepath = split(filepath, ".mat");
filename2save =  split_filepath(1)+ "_fitAug22_depthweight_1000R.mat";

%% Set up fitting settings structure
fitS.Lin2014AgeFilter = [0 4000];
fitS.weighting = "depth";    %Can be "depth", "age", or "none"
fitS.minMeanSR = 8;
fitS.minDepth = 1000;
fitS.maxLat = 40;
fitS.chi2binN = 5;
fitS.dispChi2 = false;
fitS.chi2MinCountNum = 5;   %Minimum allowed counts in a bin (chisquared)
fitS.chi2MaxCountNum = 50;  %Maximum allowed counts in a bin (chisquared)
fitS.mlnReps = 3;               %For weighting with depth, best is 3
fitS.invXbinEdges = 0:0.1:15;
fitS.enforceBinSizeLimits = true;
fitS.BMode.weightDP = 3;
fitS.BMode.weightInflator = 1;
fitS.OneRun.numruns = 1000;
if fitS.weighting == "depth"
    fitS.OneRun.weightRepDP = 3;       %For weighting with depth, best is 3
fitS.OneRun.weightRepInflator = 4; %For weighting with depth, best is 4
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
d.S1.fitS = fitS;

fitS.useParallelComp = true;
%% Fit Mix Log Norm to BIGMACS data for comparison reasons
%Load BIGMACS files
BMlognorm = readtable("../BIGMACSdata/lognormal.txt");                  % (use x values currently used in BIGMACS)
BM.nSR.x = BMlognorm.Var1';
BM.nSR.px = BMlognorm.Var2;
BM.hist = readmatrix("../BIGMACSdata/Lin2014_sedrateratio_cm_wo_NaN.txt");
BM.hist = BM.hist(:,4);
BM.TM = readmatrix("../BIGMACSdata/transition_parameter.txt");
[MLN_BIGMACS, ~, BM.gmfit] = fitMixLogNorm(BM.hist, BM.nSR.x, 2, 3);
[BM.lnSR.x, BM.lnSR.px] = px_to_pfx(BM.nSR.x, BM.nSR.px, @log);

%% Analyse new data
%Choose subset of cores
d.S1.chooseLog = d.dataT.meanSR > fitS.minMeanSR & d.dataT.depths > fitS.minDepth & abs(d.dataT.lats) < fitS.maxLat;
%Find number of cores
d.S1.numCores = sum(d.S1.chooseLog);

%% Fit MLN and inverse Gamma distribution to my Bchron Mode data
%Set up specific fitting settings
fitS.Lin2014AgeFiltering = true;

%Create x vector that spans very wide range of nSR
d.logx = -6:0.01:6;
d.x = exp(d.logx);

[d.S1.BMode] = ARfitdists(d.dataT.bchronMode, d.x, d.S1.chooseLog, fitS.BMode.weightDP, fitS.BMode.weightInflator, 1, fitS);
%% Fit dists to all BSamp samples & get chi-squared goodness of fit
countDivisor = 1000;
[d.S1.BChAR] = ARfitdists(d.dataT.bchronProb, d.x, d.S1.chooseLog, 3, 1, 1000, fitS);

%% Get Chi-Squared goodness of fits
% Do chi2gof testing on Bchron Mode data

%Make a column cell-vector that holds all pdfs to test
pdfs = {d.S1.BMode.LN.lnSR; d.S1.BMode.MLN.lnSR; d.S1.BMode.Gam.lnSR; d.S1.BMode.invGam.lnSR};
%pdfs = {d.S1.BMode.MLN.lnSR; d.S1.BMode.invGam.lnSR};
fitS.dispChi2 = true;
% h = NaN(1,1);
% p = NaN(1,1);
[h, p, chiStat] = chi2_dataVStwopdfVECs(log(d.S1.BMode.weightedC), d.S1.BMode.numCpairs, 20, pdfs, fitS);

%store the chiStats in the data structures
for i = 1:size(pdfs,1)
    chiStat{i}.h = h(i);
    chiStat{i}.p = p(i);
end
d.S1.BMode.LN.chiStats = chiStat{1};
d.S1.BMode.MLN.chiStats = chiStat{2};
d.S1.BMode.Gam.chiStats = chiStat{3};
d.S1.BMode.invGam.chiStats = chiStat{4};
%% Fit MLNs and InvGams to individual Bchron Samplings
% %Fit many individual runs
fitS.Lin2014AgeFiltering = true;
numruns = fitS.OneRun.numruns;
fitS.dispChi2 = false;
rng(2)
[d.S1.BChIR] = SRun_MLNandInvGam(d.dataT.bchronProb, d.S1.chooseLog, numruns, d.x, 2, fitS);

%% Fit dists to individual Newall Samplings
fitS.Lin2014AgeFiltering = false;
fitS.dispChi2 = false;
rng(2)
[d.S1.New0IR] = SRun_MLNandInvGam(d.dataT.nSRcounts, d.S1.chooseLog, numruns, d.x, 2, fitS);
[d.S1.New500IR] = SRun_MLNandInvGam(d.dataT.nSRcounts500, d.S1.chooseLog, numruns, d.x, 2, fitS);
[d.S1.New1000IR] = SRun_MLNandInvGam(d.dataT.nSRcounts1000, d.S1.chooseLog, numruns, d.x, 2, fitS);
[d.S1.New1500IR] = SRun_MLNandInvGam(d.dataT.nSRcounts1500, d.S1.chooseLog, numruns, d.x, 2, fitS);

%% Fit dists to collected Newall Samplings
%weight replicator = 1e-1 keeps the mean and variance the same as NewIR,
%but is faster than weightReplicator = 1;
[d.S1.New0AR] = ARfitdists(d.dataT.nSRcounts, d.x, d.S1.chooseLog, 3,1e-1, countDivisor, fitS); 
[d.S1.New500AR] = ARfitdists(d.dataT.nSRcounts500, d.x, d.S1.chooseLog, 3, 1e-1,  countDivisor, fitS);
[d.S1.New1000AR] = ARfitdists(d.dataT.nSRcounts1000, d.x, d.S1.chooseLog, 3, 1e-1, countDivisor, fitS);
[d.S1.New1500AR] = ARfitdists(d.dataT.nSRcounts1500, d.x, d.S1.chooseLog, 3, 1e-1, countDivisor, fitS);

%% Calculate histogram counts in log Space for each individual run
bw           = 0.1; %bin width
logBinEdges  = -5:bw:5;          %lognSR
d.S1.BChIR.lnSRHistCounts = NaN(numruns, length(logBinEdges)-1);
d.S1.New0IR.lnSRHistCounts = NaN(numruns, length(logBinEdges)-1);
d.S1.New500IR.lnSRHistCounts = NaN(numruns, length(logBinEdges)-1);
d.S1.New1000IR.lnSRHistCounts = NaN(numruns, length(logBinEdges)-1);
d.S1.New1500IR.lnSRHistCounts = NaN(numruns, length(logBinEdges)-1);

for i = 1:numruns
    d.S1.BChIR.lnSRHistCounts(i,:) = histcounts(log(d.S1.BChIR.weightedC{i}), 'BinEdges', logBinEdges)./d.S1.fitS.OneRun.weightRepInflator;
    d.S1.New0IR.lnSRHistCounts(i,:) = histcounts(log(d.S1.New0IR.weightedC{i}), 'BinEdges', logBinEdges)./d.S1.fitS.OneRun.weightRepInflator;
    d.S1.New500IR.lnSRHistCounts(i,:) = histcounts(log(d.S1.New500IR.weightedC{i}), 'BinEdges', logBinEdges)./d.S1.fitS.OneRun.weightRepInflator;
    d.S1.New1000IR.lnSRHistCounts(i,:) = histcounts(log(d.S1.New1000IR.weightedC{i}), 'BinEdges', logBinEdges)./d.S1.fitS.OneRun.weightRepInflator;
    d.S1.New1500IR.lnSRHistCounts(i,:) = histcounts(log(d.S1.New1500IR.weightedC{i}), 'BinEdges', logBinEdges)./d.S1.fitS.OneRun.weightRepInflator;
end

%% Chi2gof testing of Bchron Sampling data to my BchronMode fits

% Test the chi2gof of each run's weighted counts to the Bchron Mode fit
h = NaN(numruns,1);
p = NaN(numruns,1);
fitS.dispChi2 = false;
chiStat1RunT = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
for i = 1:numruns
    if ~isempty(d.S1.BChIR.weightedC{i})
        if i <4
            fitS.dispChi2 = true;
        else
            fitS.dispChi2 = false;
        end
        [h(i), p(i), chiStat1Run] = chi2gof_vsMLN(d.S1.BMode.gmfit, log(d.S1.BChIR.weightedC{i}), d.S1.BChIR.numCpairs(i), fitS);
        chiStat1RunT = [chiStat1RunT; struct2table(chiStat1Run)]; %#ok<AGROW>
    else
        emptyRow = table('Size', [1,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
        chiStat1RunT = [chiStat1RunT; emptyRow]; %#ok<AGROW>
    end
end

    %Save the results in a table and add to structure
chiStat1RunT = addvars(chiStat1RunT, h, p, 'Before', "chi2stat");
d.S1.BChIR.chiStatTvsBMode = chiStat1RunT;

%% Chi2gof testing of BMode data to BIGMACS fits
% Test the chi2gof of each fit to BIGMACS fit
fitS.dispChi2 = true;

[h, p, chiStat1Run_BM] = chi2gof_vsMLN(BM.gmfit, log(d.S1.BMode.weightedC), d.S1.BMode.numCpairs, fitS);
chiStat1RunT_BM = struct2table(chiStat1Run_BM);

%Save results to table and add to structure
chiStat1RunT_BM = addvars(chiStat1RunT_BM, h, p, 'Before', "chi2stat");
d.S1.BMode.chiStatTvsBM = chiStat1RunT_BM;

%% Chi2gof testing of Bchron Sampling data to BIGMACS fits
% Test the chi2gof of each fit to BIGMACS fit
h = NaN(numruns,1);
p = NaN(numruns,1);
fitS.dispChi2 = false;
chiStat1RunT_BM = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
for i = 1:numruns
    if ~isempty(d.S1.BChIR.weightedC{i})
        [h(i), p(i), chiStat1Run_BM] = chi2gof_vsMLN(BM.gmfit, log(d.S1.BChIR.weightedC{i}), d.S1.BChIR.numCpairs(i), fitS);
        chiStat1RunT_BM = [chiStat1RunT_BM; struct2table(chiStat1Run_BM)]; %#ok<AGROW>
    else
        emptyRow = table('Size', [1,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
        chiStat1RunT_BM = [chiStat1RunT_BM; emptyRow]; %#ok<AGROW>
    end
end

%Save results to table and add to structure
chiStat1RunT_BM = addvars(chiStat1RunT_BM, h, p, 'Before', "chi2stat");
d.S1.BChIR.chiStatTvsBM = chiStat1RunT_BM;

%% Chi2gof testing of RSR0 Sampling data to BIGMACS fits
% Test the chi2gof of each fit to BIGMACS fit
h = NaN(numruns,1);
p = NaN(numruns,1);
fitS.dispChi2 = false;
chiStat1RunT_BM = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
for i = 1:numruns
    if ~isempty(d.S1.New0IR.weightedC{i})
        [h(i), p(i), chiStat1Run_BM] = chi2gof_vsMLN(BM.gmfit, log(d.S1.New0IR.weightedC{i}), d.S1.New0IR.numCpairs(i), fitS);
        chiStat1RunT_BM = [chiStat1RunT_BM; struct2table(chiStat1Run_BM)]; %#ok<AGROW>
    else
        emptyRow = table('Size', [1,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
        chiStat1RunT_BM = [chiStat1RunT_BM; emptyRow]; %#ok<AGROW>
    end
end

%Save results to table and add to structure
chiStat1RunT_BM = addvars(chiStat1RunT_BM, h, p, 'Before', "chi2stat");
d.S1.New0IR.chiStatTvsBM = chiStat1RunT_BM;

%% Chi2gof testing of New500 Sampling data to BIGMACS fits
% Test the chi2gof of each fit to BIGMACS fit
h = NaN(numruns,1);
p = NaN(numruns,1);
fitS.dispChi2 = false;
chiStat1RunT_BM = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
for i = 1:numruns
    if ~isempty(d.S1.New500IR.weightedC{i})
        [h(i), p(i), chiStat1Run_BM] = chi2gof_vsMLN(BM.gmfit, log(d.S1.New500IR.weightedC{i}), d.S1.New500IR.numCpairs(i), fitS);
        chiStat1RunT_BM = [chiStat1RunT_BM; struct2table(chiStat1Run_BM)]; %#ok<AGROW>
    else
        emptyRow = table('Size', [1,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
        chiStat1RunT_BM = [chiStat1RunT_BM; emptyRow]; %#ok<AGROW>
    end
end

%Save results to table and add to structure
chiStat1RunT_BM = addvars(chiStat1RunT_BM, h, p, 'Before', "chi2stat");
d.S1.New500IR.chiStatTvsBM = chiStat1RunT_BM;

%% Chi2gof testing of New1000 Sampling data to BIGMACS fits
% Test the chi2gof of each fit to BIGMACS fit
h = NaN(numruns,1);
p = NaN(numruns,1);
fitS.dispChi2 = false;
chiStat1RunT_BM = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
for i = 1:numruns
    if ~isempty(d.S1.New1000IR.weightedC{i})
        [h(i), p(i), chiStat1Run_BM] = chi2gof_vsMLN(BM.gmfit, log(d.S1.New1000IR.weightedC{i}), d.S1.New1000IR.numCpairs(i), fitS);
        chiStat1RunT_BM = [chiStat1RunT_BM; struct2table(chiStat1Run_BM)]; %#ok<AGROW>
    else
        emptyRow = table('Size', [1,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
        chiStat1RunT_BM = [chiStat1RunT_BM; emptyRow]; %#ok<AGROW>
    end
end

%Save results to table and add to structure
chiStat1RunT_BM = addvars(chiStat1RunT_BM, h, p, 'Before', "chi2stat");
d.S1.New1000IR.chiStatTvsBM = chiStat1RunT_BM;

%% Chi2gof testing of New1500 Sampling data to BIGMACS fits
% Test the chi2gof of each fit to BIGMACS fit
h = NaN(numruns,1);
p = NaN(numruns,1);
fitS.dispChi2 = false;
chiStat1RunT_BM = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
for i = 1:numruns
    if ~isempty(d.S1.New1500IR.weightedC{i})
        [h(i), p(i), chiStat1Run_BM] = chi2gof_vsMLN(BM.gmfit, log(d.S1.New1500IR.weightedC{i}), d.S1.New1500IR.numCpairs(i), fitS);
        chiStat1RunT_BM = [chiStat1RunT_BM; struct2table(chiStat1Run_BM)]; %#ok<AGROW>
    else
        emptyRow = table('Size', [1,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
        chiStat1RunT_BM = [chiStat1RunT_BM; emptyRow]; %#ok<AGROW>
    end
end

%Save results to table and add to structure
chiStat1RunT_BM = addvars(chiStat1RunT_BM, h, p, 'Before', "chi2stat");
d.S1.New1500IR.chiStatTvsBM = chiStat1RunT_BM;

%% Save data as file
if isfile(filename2save)
    disp("Warning, not saving file because it already exists. Please rename file")
else
    save(filename2save, "d", "-v7.3")
end