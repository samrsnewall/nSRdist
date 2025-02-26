%% Analyse a given set of results
%Add important paths
addpath('../Functions')
addpath('../Results')

%% Load the results
filepath = "../Results/dataT_PFandLin_R200M20_Feb19_nSRmeanTrue";
load(filepath)
d.dataT = dataT;
d.S = S;
d.label = "AllCores, R200M20, nSRmTrue";
split_filepath = split(filepath, ".mat");
filename2save =  split_filepath(1)+ "_fitFeb19.mat";

%% Set up fitting settings structure
fitS.Lin2014AgeFiltering = 1;
fitS.weighting = "depth";
fitS.chi2binN = 5;
fitS.dispChi2 = false;
fitS.mlnReps = 5;
fitS.invXbinEdges = 0:0.1:15;
fitS.enforceBinSizeLimits = true;
fitS.OneRun.weightRepDP = 3;
fitS.OneRun.weightRepInflator = 4;
fitS.OneRun.MLNReps = 3;
d.S1.fitS = fitS;

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
d.S1.chooseLog = d.dataT.meanSR > 8;
%Find number of cores
d.S1.numCores = sum(d.S1.chooseLog);

%% Fit MLN and inverse Gamma distribution to my Bchron Mode data
%Create x vector that spans very wide range of nSR
logx = -6:0.01:6;
d.x = exp(logx);

disp([d.label;  "BchronMode vs best fit mln"])
[d.S1.BMode.mixLog, d.S1.BMode.Hist,~,~,~,...
    d.S1.BMode.gmfit, d.S1.BMode.nc, d.S1.BMode.h, d.S1.BMode.p,d.S1.BMode.chiStat]...
    = plotSRandResHistograms(dataT.bchronMode,...
    d.x, d.S1.chooseLog, 3, 1, 2, 0, "", 0, fitS);
if fitS.dispChi2
    gcf; title([d.label; "chi2gof: Bmode Data vs Best Fit MLN"])
end

d.S1.BMode.MLN.nSR.x = d.S1.BMode.mixLog(:,1);
d.S1.BMode.MLN.nSR.px = d.S1.BMode.mixLog(:,2);

%Fit Inverse gamma distribution to my Bchron Mode data
[d.S1.BMode.invGam.nSR.x, d.S1.BMode.invGam.nSR.px] = fitGamma2invSR(d.dataT.bchronMode, d.S1.chooseLog, fitS);

% Do chi2gof testing on Bchron Mode data
%Need to change inverse gamma pdf to log(nSR) space
f_log = @(x) log(x);
[d.S1.BMode.MLN.lnSR.x, d.S1.BMode.MLN.lnSR.px] = px_to_pfx(d.S1.BMode.MLN.nSR.x, d.S1.BMode.MLN.nSR.px, f_log);
d.S1.BMode.pdfVEC1.x = d.S1.BMode.MLN.lnSR.x;
d.S1.BMode.pdfVEC1.px = d.S1.BMode.MLN.lnSR.px;
d.S1.BMode.pdfVEC1.numParams = 6;

[d.S1.BMode.invGam.lnSR.x, d.S1.BMode.invGam.lnSR.px] = px_to_pfx(d.S1.BMode.invGam.nSR.x, d.S1.BMode.invGam.nSR.px, f_log);
d.S1.BMode.pdfVEC2.x = d.S1.BMode.invGam.lnSR.x;
d.S1.BMode.pdfVEC2.px = d.S1.BMode.invGam.lnSR.px;
d.S1.BMode.pdfVEC2.numParams = 2;

%%
fitS.dispChi2 = true;
[h1, p1, chiStat1, h2, p2, chiStat2] = chi2_dataVStwopdfVECs(log(d.S1.BMode.Hist), d.S1.BMode.nc, 20, d.S1.BMode.pdfVEC1, d.S1.BMode.pdfVEC2, fitS);
chiStat1.h = h1;
chiStat1.p = p1;
chiStat2.h = h2;
chiStat2.p = p2;
d.S1.BMode.MLN.chiStat = chiStat1;
d.S1.BMode.invGam.chiStat = chiStat2;

%% Fit MLNs and InvGams to individual Bchron Samplings
% %Fit many individual runs
numruns = 100;
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

%% Chi2gof testing of Bchron Sampling data to other fits

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

if isfile(filename2save)
    disp("Warning, not saving file because it already exists, rename file")
else
    save(filename2save, "d")
end