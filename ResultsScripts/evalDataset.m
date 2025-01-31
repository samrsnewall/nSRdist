%% Analyse a given set of results
%Add important paths
addpath('../Functions')

%Load the results
load("../Results/dataT_LinandPF_R200M20_Dec12.mat")
d.dataT = dataT;
d.S = S;

%Set up fitting settings structure
fitS.Lin2014AgeFiltering = 1;
fitS.weighting = "depth";
fitS.chi2binN = 10;
fitS.dispChi2 = false;
fitS.mln1RunReps = 1;
fitS.mlnReps = 5;
fitS.enforceBinSizeLimits = true;

%% Fit Mix Log Norm to BIGMACS data for comparison reasons
%Load BIGMACS files
lognorm_BIGMACS = readtable("../BIGMACSdata/lognormal.txt");                  % (use x values currently used in BIGMACS)
x = lognorm_BIGMACS.Var1';
BIGMACShist = readmatrix("../BIGMACSdata/Lin2014_sedrateratio_cm_wo_NaN.txt");
BIGMACShist = BIGMACShist(:,4);
TM_BIGMACS = readmatrix("../BIGMACSdata/transition_parameter.txt");
[MLN_BIGMACS, ~, gmfitBM] = fitMixLogNorm(BIGMACShist, x, 2, 0, fitS.mlnReps);

%% Analyse new data
%Find number of cores
sizeCores = numel(d.dataT.lats);
numCores = sizeCores;
%Find high SR subset
highSRLog = d.dataT.meanSR > 8;

%% Output Summary Figures
outputMetadataAndSummaryFigures(highSRLog,d.dataT)

%% Fit MLN distribution to my Bchron Mode data
disp("LinMethod (mode) data vs best fit mln")
[d.mixLogBMode, d.BModeHist,~,~,~,d.gmfitBmode, d.ncBmode, d.h, d.p,d.chiStat] = plotSRandResHistograms(dataT.bchronMode, x, highSRLog, 3, 1, 2, 0, "", 1, fitS);
if fitS.dispChi2
    gcf; title("chi2gof: NoLinCoreBmode Data vs Best Fit MLN")
end

%Plot new histogram and MLN alongside BIGMACS MLN
figure;
hold on
yyaxis("left")
histogram(d.BModeHist, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts')
ylabel("Counts")
yyaxis("right")
plot(x, lognorm_BIGMACS.Var2, '-k', 'DisplayName', "BIGMACS", 'LineWidth', 1)
plot(x,d.mixLogBMode(:,2), '-r', 'DisplayName', "Replicated MLN", 'LineWidth', 1)
xlim([0 6])
xlabel("nSR")
ylabel("PDF")
legend()
title("BIGMACS replicated with independent cores")
set(gcf, 'Position', [20, 400, 500, 200])

%See how my distribution performs with chi2gof on BIGMACS gmfit
disp("LinMethod data (mode) vs BIGMACS MLN")
[d.hBM, d.pBM, d.chiStatBM] = chi2gof_vsMLN(gmfitBM, log(d.BModeHist), d.ncBmode, fitS);
if fitS.dispChi2
    gcf; title("chi2gof: NoLinCoreBmode Data vs BIGMACS fit")
end

%% Fit MLNs to individual Bchron Samplings
%Fit many individual runs
numruns = 100;
fitS.dispChi2 = false;
rng(2)
[d.MLN1R.pdfs, d.MLN1R.c95up,...
    d.MLN1R.c95down, d.MLN1R.mus,...
    d.MLN1R.sigmas, d.MLN1R.outputS]...
    = SingleRunLogNorms(d.dataT.bchronProb,true(sizeCores,1),...
    numruns, x, 2, 3, 4, 0, fitS);

%% Fit Inverse Gammas to individual Bchron Samplings
[d.InvGam1R.invSRpdfs, ~,~,d.InvGam1R.phats, d.InvGam1R.nSRpdfs] = SingleRunGammas(d.dataT.bchronProb,true(sizeCores,1),...
    numruns, x, 3, 4, 0, fitS);

%% Plot a histogram of all the weighted counts

%First, plot histogram of BchronMode data
figure;
subplot(3,1,1)
hBchMode = histogram(d.BModeHist, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts');
xlim([0 6])
title("BchronMode")

%Using BchronMode binEdges, plot histogram of all sampling data combined
WCall = [];
oneRunBinCounts = NaN(numruns, length(hBchMode.BinEdges)-1);
for i = 1:numruns
    WCall = [WCall, d.MLN1R.outputS.weightedC{i}];
    oneRunBinCounts(i,:) = histcounts(d.MLN1R.outputS.weightedC{i}, hBchMode.BinEdges);
end

subplot(3,1,2)
histogram(WCall, 'BinEdges', hBchMode.BinEdges, 'FaceColor', 'k', 'FaceAlpha', 0.4)
xlim([0 6])
title("Bchron Samples Combined")

%Now plot each individual histogram with a transparent face color
subplot(3,1,3)
for i = 1:numruns
    hold on
    histogram('BinCounts',oneRunBinCounts(i,:), 'BinEdges', hBchMode.BinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none')
end
xlim([0 6])
xlabel("nSR")
title("Bchron Samples, Transparent Histograms")

%% Plot histograms in log nSR space
figure;
subplot(3,1,1)
hBchModeLog = histogram(log(d.BModeHist), 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts');
xlim([-4 4])
title("BchronMode")

oneRunBinCounts_lognSR = NaN(numruns, length(hBchModeLog.BinEdges)-1);
for i = 1:numruns
    oneRunBinCounts_lognSR(i,:) = histcounts(log(d.MLN1R.outputS.weightedC{i}), hBchModeLog.BinEdges);
end
subplot(3,1,2)
histogram(log(WCall), 'BinEdges', hBchModeLog.BinEdges, 'FaceColor', 'k', 'FaceAlpha', 0.4)
xlim([-4 4])
title("Bchron Samples Combined")
subplot(3,1,3)
for i = 1:numruns
    hold on
    histogram('BinCounts',oneRunBinCounts_lognSR(i,:), 'BinEdges', hBchModeLog.BinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none')
end
xlim([-4 4])
xlabel("log(nSR)")
title("Bchron Samples, Transparent Histograms")

%% Do chi2gof testing
% Test the chi2gof of each run's fit to it's Bchron Mode fit
h1R = NaN(numruns,1);
p1R = NaN(numruns,1);
fitS.dispChi2 = false;
chiStat1RunT = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
for i = 1:numruns
    if ~isempty(d.MLN1R.outputS.weightedC{i})
    [h1R(i), p1R(i), chiStat1Run] = chi2gof_vsMLN(d.gmfitBmode, log(d.MLN1R.outputS.weightedC{i}), d.MLN1R.outputS.numCpairs(i), fitS);
    chiStat1RunT = [chiStat1RunT; struct2table(chiStat1Run)]; %#ok<AGROW>
    end
end

%Save the results in a table and add to structure
chiStat1RunT = addvars(chiStat1RunT, h1R, p1R, 'Before', "chi2stat");
d.MLN1R.chiStat1RunT = chiStat1RunT;

% Test the chi2gof of each fit to BIGMACS fit
h1R = NaN(numruns,1);
p1R = NaN(numruns,1);
fitS.dispChi2 = false;
chiStat1RunT_BM = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
for i = 1:numruns
    [h1R(i), p1R(i), chiStat1Run_BM] = chi2gof_vsMLN(gmfitBM, log(d.MLN1R.outputS.weightedC{i}), d.MLN1R.outputS.numCpairs(i), fitS);
    chiStat1RunT_BM = [chiStat1RunT_BM; struct2table(chiStat1Run_BM)]; %#ok<AGROW>
end

%Save results to table and add to structure
chiStat1RunT_BM = addvars(chiStat1RunT_BM, h1R, p1R, 'Before', "chi2stat");
d.MLN1R.chiStat1RunT_BM = chiStat1RunT_BM;

%Plot a figure of all the fits
figure;
hold on
plot(x, d.MLN1R.pdfs, 'Color', [0,0,0,0.1], 'HandleVisibility', 'off')
plot(x, d.mixLogBMode(:,2), '-r', 'DisplayName', "BchronMode: NewCores,Marine20", 'LineWidth', 1)
plot(NaN, NaN,'LineStyle', "none", 'DisplayName', "Random Runs that reject H0 = " + num2str(mean(d.MLN1R.chiStat1RunT.h1R).*100, 3) + "%")
plot(x, MLN_BIGMACS(:,2), '-k', 'DisplayName', "BchronMode: BIGMACS", 'LineWidth', 1)
plot(NaN, NaN, 'LineStyle', "none", 'DisplayName', "Random Runs that reject H0 = " + num2str(mean(d.MLN1R.chiStat1RunT_BM.h1R).*100, 3) + "%")
xlim([0 6])
%ylim(commonYLim)
legend()
xlabel("nSR")
ylabel("PDF")
title("Independent Cores; Marine20; R = 0±200")

