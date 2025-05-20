%% Analyse a given set of results
%Add important paths
addpath('../Functions')
addpath('../Results')

%% Load the results
filepath = "../Results/dataT_RLGtrueB_R200M20_May14.mat";
load(filepath)
d.dataT = dataT;
d.S = S;
d.label = "AllCores, R200M20";
split_filepath = split(filepath, ".mat");
filename2save =  split_filepath(1)+ "_fitMay14_depthweight.mat";

%% Set up fitting settings structure
fitS.Lin2014AgeFilter = [500 4000];
fitS.weighting = "none";    %Can be "depth", "age", or "none"
fitS.chi2binN = 5;
fitS.dispChi2 = false;
fitS.chi2MinCountNum = 0;
fitS.mlnReps = 3;               %For weighting with depth, best is 3
fitS.invXbinEdges = 0:0.1:15;
fitS.enforceBinSizeLimits = true;
fitS.OneRun.numruns = 400;
fitS.OneRun.weightRepDP = 3;       %For weighting with depth, best is 3
fitS.OneRun.weightRepInflator = 4; %For weighting with depth, best is 4
fitS.OneRun.MLNReps = 3;           %For weighting with depth, best is 3
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
d.S1.chooseLog = d.dataT.meanSR > 8 & d.dataT.depths > 1000 & abs(d.dataT.lats) < 40;
%Find number of cores
d.S1.numCores = sum(d.S1.chooseLog);

%% Fit MLN and inverse Gamma distribution to my Bchron Mode data

%Set up specific fitting settings
fitS.Lin2014AgeFiltering = true;

%Create x vector that spans very wide range of nSR
d.logx = -6:0.01:6;
d.x = exp(d.logx);

%Fit Mix Log Normal
disp([d.label;  "BchronMode vs best fit mln"])
[d.S1.BMode.mixLog, d.S1.BMode.weightedC,d.S1.BMode.agediffs,~,~,...
    d.S1.BMode.gmfit, d.S1.BMode.numCpairs, d.S1.BMode.sedLength, ...
    d.S1.BMode.sedTimeSpan, ~, ~, ~]...
    = plotSRandResHistograms(dataT.bchronMode,...
    d.x, d.S1.chooseLog, 3, 1, 2, 0, "", 0, fitS);
if fitS.dispChi2
    gcf; title([d.label; "chi2gof: Bmode Data vs Best Fit MLN"])
end

d.S1.BMode.MLN.nSR.x = d.S1.BMode.mixLog(:,1);
d.S1.BMode.MLN.nSR.px = d.S1.BMode.mixLog(:,2);

%Fit Inverse gamma distribution to my Bchron Mode data
[d.S1.BMode.invGam.nSR.x, d.S1.BMode.invGam.nSR.px] = fitGamma2invSR(d.dataT.bchronMode, d.S1.chooseLog, fitS);

%Fit gamma distribution to my Bchron Mode data
[d.S1.BMode.Gam.nSR.x, d.S1.BMode.Gam.nSR.px] = fitGamma2nSR(d.dataT.bchronMode, d.S1.chooseLog, fitS);

%Fit LogNormal distribution to my Bchron Mode data
[d.S1.BMode.LN.nSR.x, d.S1.BMode.LN.nSR.px] = fitLogNorm2nSR(d.dataT.bchronMode, d.S1.chooseLog, fitS);


%% Get Chi-Squared goodness of fits
% Do chi2gof testing on Bchron Mode data
%Need to change pdfs to log(nSR) space
f_log = @(x) log(x);
[d.S1.BMode.MLN.lnSR.x, d.S1.BMode.MLN.lnSR.px] = px_to_pfx(d.S1.BMode.MLN.nSR.x, d.S1.BMode.MLN.nSR.px, f_log);
d.S1.BMode.MLN.lnSR.numParams = 6;
d.S1.BMode.MLN.lnSR.pdfName = "2 Component Mix Log Normal";

[d.S1.BMode.invGam.lnSR.x, d.S1.BMode.invGam.lnSR.px] = px_to_pfx(d.S1.BMode.invGam.nSR.x, d.S1.BMode.invGam.nSR.px, f_log);
d.S1.BMode.invGam.lnSR.numParams = 2;
d.S1.BMode.invGam.lnSR.pdfName = "Inverse Gamma";

[d.S1.BMode.Gam.lnSR.x, d.S1.BMode.Gam.lnSR.px] = px_to_pfx(d.S1.BMode.Gam.nSR.x, d.S1.BMode.Gam.nSR.px, f_log);
d.S1.BMode.Gam.lnSR.numParams = 2;
d.S1.BMode.Gam.lnSR.pdfName = "Gamma";

[d.S1.BMode.LN.lnSR.x, d.S1.BMode.LN.lnSR.px] = px_to_pfx(d.S1.BMode.LN.nSR.x, d.S1.BMode.LN.nSR.px, f_log);
d.S1.BMode.LN.lnSR.numParams = 2;
d.S1.BMode.LN.lnSR.pdfName = "LogNorm";

%Make a column cell-vector that holds all pdfs to test
pdfs = {d.S1.BMode.LN.lnSR; d.S1.BMode.MLN.lnSR; d.S1.BMode.Gam.lnSR; d.S1.BMode.invGam.lnSR};
pdfs = {d.S1.BMode.MLN.lnSR; d.S1.BMode.invGam.lnSR};
fitS.dispChi2 = true;
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

%% Fit MLNs and InvGams to individual Newall Samplings with 0yr restriction
fitS.Lin2014AgeFiltering = false;
fitS.dispChi2 = false;
rng(2)
[d.S1.New0IR] = SRun_MLNandInvGam(d.dataT.nSRcounts, d.S1.chooseLog, numruns, d.x, 2, fitS);

%% Fit MLNs and InvGams to individual Newall Samplings with 500yr restriction
fitS.Lin2014AgeFiltering = false;
fitS.dispChi2 = false;
rng(2)
[d.S1.New500IR] = SRun_MLNandInvGam(d.dataT.nSRcounts500, d.S1.chooseLog, numruns, d.x, 2, fitS);

%% Fit MLNs and InvGams to individual Newall Samplings with 1000yr restriction
fitS.Lin2014AgeFiltering = false;
fitS.dispChi2 = false;
rng(2)
[d.S1.New1000IR] = SRun_MLNandInvGam(d.dataT.nSRcounts1000, d.S1.chooseLog, numruns, d.x, 2, fitS);

%% Fit MLNs and InvGams to individual Newall Samplings with 1500yr restriction
fitS.Lin2014AgeFiltering = false;
fitS.dispChi2 = false;
rng(2)
[d.S1.New1500IR] = SRun_MLNandInvGam(d.dataT.nSRcounts1500, d.S1.chooseLog, numruns, d.x, 2, fitS);

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
h1R = NaN(numruns,1);
p1R = NaN(numruns,1);
fitS.dispChi2 = false;
chiStat1RunT = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
for i = 1:numruns
    if ~isempty(d.S1.BChIR.weightedC{i})
        [h1R(i), p1R(i), chiStat1Run] = chi2gof_vsMLN(d.S1.BMode.gmfit, log(d.S1.BChIR.weightedC{i}), d.S1.BChIR.numCpairs(i), fitS);
        chiStat1RunT = [chiStat1RunT; struct2table(chiStat1Run)]; %#ok<AGROW>
    else
        emptyRow = table('Size', [1,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
        chiStat1RunT = [chiStat1RunT; emptyRow]; %#ok<AGROW>
    end
end

    %Save the results in a table and add to structure
chiStat1RunT = addvars(chiStat1RunT, h1R, p1R, 'Before', "chi2stat");
d.S1.BChIR.chiStatTvsBMode = chiStat1RunT;

%% Chi2gof testing of Bchron Sampling data to BIGMACS fits
% Test the chi2gof of each fit to BIGMACS fit
h1R = NaN(numruns,1);
p1R = NaN(numruns,1);
fitS.dispChi2 = false;
chiStat1RunT_BM = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
for i = 1:numruns
    if ~isempty(d.S1.BChIR.weightedC{i})
        [h1R(i), p1R(i), chiStat1Run_BM] = chi2gof_vsMLN(BM.gmfit, log(d.S1.BChIR.weightedC{i}), d.S1.BChIR.numCpairs(i), fitS);
        chiStat1RunT_BM = [chiStat1RunT_BM; struct2table(chiStat1Run_BM)]; %#ok<AGROW>
    else
        emptyRow = table('Size', [1,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
        chiStat1RunT_BM = [chiStat1RunT_BM; emptyRow]; %#ok<AGROW>
    end
end

%Save results to table and add to structure
chiStat1RunT_BM = addvars(chiStat1RunT_BM, h1R, p1R, 'Before', "chi2stat");
d.S1.BChIR.chiStatTvsBM = chiStat1RunT_BM;

%% Chi2gof testing of RSR0 Sampling data to BIGMACS fits
% Test the chi2gof of each fit to BIGMACS fit
h1R = NaN(numruns,1);
p1R = NaN(numruns,1);
fitS.dispChi2 = false;
chiStat1RunT_BM = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
for i = 1:numruns
    if ~isempty(d.S1.New0IR.weightedC{i})
        [h1R(i), p1R(i), chiStat1Run_BM] = chi2gof_vsMLN(BM.gmfit, log(d.S1.New0IR.weightedC{i}), d.S1.New0IR.numCpairs(i), fitS);
        chiStat1RunT_BM = [chiStat1RunT_BM; struct2table(chiStat1Run_BM)]; %#ok<AGROW>
    else
        emptyRow = table('Size', [1,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
        chiStat1RunT_BM = [chiStat1RunT_BM; emptyRow]; %#ok<AGROW>
    end
end

%Save results to table and add to structure
chiStat1RunT_BM = addvars(chiStat1RunT_BM, h1R, p1R, 'Before', "chi2stat");
d.S1.New0IR.chiStatTvsBM = chiStat1RunT_BM;

%% Chi2gof testing of New500 Sampling data to BIGMACS fits
% Test the chi2gof of each fit to BIGMACS fit
h1R = NaN(numruns,1);
p1R = NaN(numruns,1);
fitS.dispChi2 = false;
chiStat1RunT_BM = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
for i = 1:numruns
    if ~isempty(d.S1.New500IR.weightedC{i})
        [h1R(i), p1R(i), chiStat1Run_BM] = chi2gof_vsMLN(BM.gmfit, log(d.S1.New500IR.weightedC{i}), d.S1.New500IR.numCpairs(i), fitS);
        chiStat1RunT_BM = [chiStat1RunT_BM; struct2table(chiStat1Run_BM)]; %#ok<AGROW>
    else
        emptyRow = table('Size', [1,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
        chiStat1RunT_BM = [chiStat1RunT_BM; emptyRow]; %#ok<AGROW>
    end
end

%Save results to table and add to structure
chiStat1RunT_BM = addvars(chiStat1RunT_BM, h1R, p1R, 'Before', "chi2stat");
d.S1.New500IR.chiStatTvsBM = chiStat1RunT_BM;

%% Chi2gof testing of New1000 Sampling data to BIGMACS fits
% Test the chi2gof of each fit to BIGMACS fit
h1R = NaN(numruns,1);
p1R = NaN(numruns,1);
fitS.dispChi2 = false;
chiStat1RunT_BM = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
for i = 1:numruns
    if ~isempty(d.S1.New1000IR.weightedC{i})
        [h1R(i), p1R(i), chiStat1Run_BM] = chi2gof_vsMLN(BM.gmfit, log(d.S1.New1000IR.weightedC{i}), d.S1.New1000IR.numCpairs(i), fitS);
        chiStat1RunT_BM = [chiStat1RunT_BM; struct2table(chiStat1Run_BM)]; %#ok<AGROW>
    else
        emptyRow = table('Size', [1,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
        chiStat1RunT_BM = [chiStat1RunT_BM; emptyRow]; %#ok<AGROW>
    end
end

%Save results to table and add to structure
chiStat1RunT_BM = addvars(chiStat1RunT_BM, h1R, p1R, 'Before', "chi2stat");
d.S1.New1000IR.chiStatTvsBM = chiStat1RunT_BM;

%% Chi2gof testing of New1500 Sampling data to BIGMACS fits
% Test the chi2gof of each fit to BIGMACS fit
h1R = NaN(numruns,1);
p1R = NaN(numruns,1);
fitS.dispChi2 = false;
chiStat1RunT_BM = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
for i = 1:numruns
    if ~isempty(d.S1.New1500IR.weightedC{i})
        [h1R(i), p1R(i), chiStat1Run_BM] = chi2gof_vsMLN(BM.gmfit, log(d.S1.New1500IR.weightedC{i}), d.S1.New1500IR.numCpairs(i), fitS);
        chiStat1RunT_BM = [chiStat1RunT_BM; struct2table(chiStat1Run_BM)]; %#ok<AGROW>
    else
        emptyRow = table('Size', [1,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
        chiStat1RunT_BM = [chiStat1RunT_BM; emptyRow]; %#ok<AGROW>
    end
end

%Save results to table and add to structure
chiStat1RunT_BM = addvars(chiStat1RunT_BM, h1R, p1R, 'Before', "chi2stat");
d.S1.New1500IR.chiStatTvsBM = chiStat1RunT_BM;

%% Save data as file
if isfile(filename2save)
    disp("Warning, not saving file because it already exists. Please rename file")
else
    save(filename2save, "d", "-v7.3")
end