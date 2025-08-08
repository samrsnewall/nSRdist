%% REPLICATING BIGMACS MLN - Using WA2022 radiocarbon cores that were NOT in Lin2014 (and are >1000m and <40deg N or S latitudinally, and have >4 radiocarbon dates)
%Add important paths
addpath('../Functions')

%Load the results
%load("../Results/dataT_LinandPF_R200M20_Dec12.mat")
load("../Results/dataT_Independent2_R200M20_Dec19.mat")
NoLinCoreRep.dataT = dataT;
NoLinCoreRep.S = S;

%Load BIGMACS files
lognorm_BIGMACS = readtable("../BIGMACSdata/lognormal.txt");                  % (use x values currently used in BIGMACS)
x = lognorm_BIGMACS.Var1';
BIGMACShist = readmatrix("../BIGMACSdata/Lin2014_sedrateratio_cm_wo_NaN.txt");
BIGMACShist = BIGMACShist(:,4);
TM_BIGMACS = readmatrix("../BIGMACSdata/transition_parameter.txt");

%Set up fitting settings structure
fitS.Lin2014AgeFiltering = 1;
fitS.Lin2014AgeFilter = [500 4000];
fitS.weighting = "depth";
fitS.chi2binN = 10;
fitS.dispChi2 = true;
fitS.mln1RunReps = 1;
fitS.mlnReps = 5;
fitS.enforceBinSizeLimits = false;

%Fit Mix Log Norm to data
[MLN_BIGMACS, ~, gmfitBM] = fitMixLogNorm(BIGMACShist, x, 2, fitS.mlnReps);

%Find number of cores
sizeCores = numel(NoLinCoreRep.dataT.lats);
numCores = sizeCores;

%Find high SR subset
highSRLog = NoLinCoreRep.dataT.meanSR > 8;

outputMetadataAndSummaryFigures(highSRLog,NoLinCoreRep.dataT)

%Fit MLN distribution to my Bchron Mode data
disp("LinMethod (mode) data vs best fit mln")
[NoLinCoreRep.mixLogBMode, NoLinCoreRep.BModeHist,~,~,~,NoLinCoreRep.gmfitBmode, NoLinCoreRep.ncBmode, NoLinCoreRep.h, NoLinCoreRep.p,NoLinCoreRep.chiStat] = plotSRandResHistograms(dataT.bchronMode, x, highSRLog, 3, 1, 2, 0, "", 1, fitS);
if fitS.dispChi2
    gcf; title("chi2gof: NoLinCoreBmode Data vs Best Fit MLN")
end

%Plot new histogram and MLN alongside BIGMACS MLN
figure;
hold on
yyaxis("left")
histogram(NoLinCoreRep.BModeHist, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts')
ylabel("Counts")
yyaxis("right")
plot(x, lognorm_BIGMACS.Var2, '-k', 'DisplayName', "BIGMACS", 'LineWidth', 1)
plot(x,NoLinCoreRep.mixLogBMode(:,2), '-r', 'DisplayName', "Replicated MLN", 'LineWidth', 1)
xlim([0 6])
xlabel("nSR")
ylabel("PDF")
legend()
title("BIGMACS replicated with independent cores")
set(gcf, 'Position', [20, 400, 500, 200])

%See how my distribution performs with chi2gof on BIGMACS gmfit
disp("LinMethod data (mode) vs BIGMACS MLN")
[NoLinCoreRep.hBM, NoLinCoreRep.pBM, NoLinCoreRep.chiStatBM] = chi2gof_vsMLN(gmfitBM, log(NoLinCoreRep.BModeHist), NoLinCoreRep.ncBmode, fitS);
if fitS.dispChi2
    gcf; title("chi2gof: NoLinCoreBmode Data vs BIGMACS fit")
end


%%%% OLD COMMENTS BELOW (When I run this again I don't see the same results as
%%%% suggested below

%The distribution here doesn't pass the chi2gof test for its own best-fit
%distribution. One core that contributes to this is probably MD06-3075,
%which has an nSR of 5.9 for 228 cm. Authors believe this is due to
%anomalous sediment oversampling during the coring process, as opposed to
%being due to a greater sedimentation rate. What happens without this core?

%%%% OLD COMMENTS ABOVE

%% REPLICATING BIGMACS MLN - Using WA2022 radiocarbon cores that were NOT in Lin2014 (and are >1000m and <40deg N or S latitudinally, and have >4 radiocarbon dates)

%Load the results again (just to ensure correct start point)
load("../Results/dataT_Independent2_R200M20_Dec19.mat")

%Remove core(s)
logMD06_3075 = ismember(dataT.cores, 'MD06-3075');
NoLinCoreRep1.dataT = dataT(~logMD06_3075, :);
NoLinCoreRep1.S = S;

%Find number of cores
sizeCores = numel(NoLinCoreRep1.dataT.lats);
numCores = sizeCores;

%Fit MLN distribution to my Bchron Mode data
disp("LinMethod (mode) data vs best fit mln")
[NoLinCoreRep1.mixLogBMode, NoLinCoreRep1.BModeHist,~,~,~,NoLinCoreRep1.gmfitBmode, NoLinCoreRep1.ncBmode, NoLinCoreRep1.h, NoLinCoreRep1.p,NoLinCoreRep1.chiStat] = plotSRandResHistograms(NoLinCoreRep1.dataT.bchronMode, x, true(sizeCores), 3, 1, 2, 0, "", 1, fitS);
if fitS.dispChi2
    gcf; title("chi2gof: NoLinCoreBmode Data vs Best Fit MLN")
end

%Plot new histogram and MLN alongside BIGMACS MLN
figure;
hold on
yyaxis("left")
histogram(NoLinCoreRep1.BModeHist, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts')
ylabel("Counts")
yyaxis("right")
plot(x, lognorm_BIGMACS.Var2, '-k', 'DisplayName', "BIGMACS", 'LineWidth', 1)
plot(x,NoLinCoreRep1.mixLogBMode(:,2), '-r', 'DisplayName', "Replicated MLN", 'LineWidth', 1)
xlim([0 6])
xlabel("nSR")
ylabel("PDF")
legend()
title("BIGMACS replicated with independent cores")
set(gcf, 'Position', [20, 400, 500, 200])

%See how my distribution performs with chi2gof on BIGMACS gmfit
disp("LinMethod data (mode) vs BIGMACS MLN")
[NoLinCoreRep1.hBM, NoLinCoreRep1.pBM, NoLinCoreRep1.chiStatBM] = chi2gof_vsMLN(gmfitBM, log(NoLinCoreRep1.BModeHist), NoLinCoreRep1.ncBmode, fitS);
if fitS.dispChi2
    gcf; title("chi2gof: NoLinCoreBmode Data vs BIGMACS fit")
end

%% Try getting chi2gof for the individual Bchron runs

%Fit many individual runs
numruns = 100;
fitS.dispChi2 = false;
rng(2)
[NoLinCoreRep.MLN1R.pdfs, NoLinCoreRep.MLN1R.c95up,...
    NoLinCoreRep.MLN1R.c95down, NoLinCoreRep.MLN1R.mus,...
    NoLinCoreRep.MLN1R.sigmas, NoLinCoreRep.MLN1R.outputS]...
    = SingleRunLogNorms(NoLinCoreRep.dataT.bchronProb,true(sizeCores,1),...
    numruns, x, 2, 3, 4, 0, fitS);

%Plot a histogram of all the weighted counts
figure;
subplot(3,1,1)
hBchMode = histogram(NoLinCoreRep.BModeHist, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts');
xlim([0 6])
title("BchronMode")

WCall = [];
oneRunBinCounts = NaN(numruns, length(hBchMode.BinEdges)-1);
for i = 1:numruns
    WCall = [WCall, NoLinCoreRep.MLN1R.outputS.weightedC{i}];
    oneRunBinCounts(i,:) = histcounts(NoLinCoreRep.MLN1R.outputS.weightedC{i}, hBchMode.BinEdges);
end

figure;
subplot(3,1,1)
hBchMode = histogram(NoLinCoreRep.BModeHist, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts');
xlim([0 6])
title("BchronMode")
subplot(3,1,2)
histogram(WCall, 'BinEdges', hBchMode.BinEdges, 'FaceColor', 'k', 'FaceAlpha', 0.4)
xlim([0 6])
title("Bchron Samples Combined")
subplot(3,1,3)
for i = 1:numruns
    hold on
    histogram('BinCounts',oneRunBinCounts(i,:), 'BinEdges', hBchMode.BinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none')
end
xlim([0 6])
xlabel("nSR")
title("Bchron Samples, Transparent Histograms")

figure;
subplot(3,1,1)
hBchModeLog = histogram(log(NoLinCoreRep.BModeHist), 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts');
xlim([-4 4])
title("BchronMode")

oneRunBinCounts_lognSR = NaN(numruns, length(hBchModeLog.BinEdges)-1);
for i = 1:numruns
    oneRunBinCounts_lognSR(i,:) = histcounts(log(NoLinCoreRep.MLN1R.outputS.weightedC{i}), hBchModeLog.BinEdges);
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

% Test the chi2gof of each run's fit to Itself
h1R = NaN(numruns,1);
p1R = NaN(numruns,1);
fitS.dispChi2 = false;
chiStat1RunT = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
for i = 1:numruns
    if ~isempty(NoLinCoreRep.MLN1R.outputS.weightedC{i})
    [h1R(i), p1R(i), chiStat1Run] = chi2gof_vsMLN(NoLinCoreRep.gmfitBmode, log(NoLinCoreRep.MLN1R.outputS.weightedC{i}), NoLinCoreRep.MLN1R.outputS.numCpairs(i), fitS);
    chiStat1RunT = [chiStat1RunT; struct2table(chiStat1Run)]; %#ok<AGROW>
    end
end

%Save the results in a table and add to structure
chiStat1RunT = addvars(chiStat1RunT, h1R, p1R, 'Before', "chi2stat");
NoLinCoreRep.MLN1R.chiStat1RunT = chiStat1RunT;

% Test the chi2gof of each fit to BIGMACS
h1R = NaN(numruns,1);
p1R = NaN(numruns,1);
fitS.dispChi2 = false;
chiStat1RunT_BM = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
for i = 1:numruns
    [h1R(i), p1R(i), chiStat1Run_BM] = chi2gof_vsMLN(gmfitBM, log(NoLinCoreRep.MLN1R.outputS.weightedC{i}), NoLinCoreRep.MLN1R.outputS.numCpairs(i), fitS);
    chiStat1RunT_BM = [chiStat1RunT_BM; struct2table(chiStat1Run_BM)]; %#ok<AGROW>
end

%Save results to table and add to structure
chiStat1RunT_BM = addvars(chiStat1RunT_BM, h1R, p1R, 'Before', "chi2stat");
NoLinCoreRep.MLN1R.chiStat1RunT_BM = chiStat1RunT_BM;

%Plot a figure of all the fits
figure;
hold on
plot(x, NoLinCoreRep.MLN1R.pdfs, 'Color', [0,0,0,0.1], 'HandleVisibility', 'off')
plot(x, NoLinCoreRep.mixLogBMode(:,2), '-r', 'DisplayName', "BchronMode: NewCores,Marine20", 'LineWidth', 1)
plot(NaN, NaN,'LineStyle', "none", 'DisplayName', "Random Runs that reject H0 = " + num2str(mean(NoLinCoreRep.MLN1R.chiStat1RunT.h1R).*100, 3) + "%")
plot(x, MLN_BIGMACS(:,2), '-k', 'DisplayName', "BchronMode: BIGMACS", 'LineWidth', 1)
plot(NaN, NaN, 'LineStyle', "none", 'DisplayName', "Random Runs that reject H0 = " + num2str(mean(NoLinCoreRep.MLN1R.chiStat1RunT_BM.h1R).*100, 3) + "%")
xlim([0 6])
%ylim(commonYLim)
legend()
xlabel("nSR")
ylabel("PDF")
title("Independent Cores; Marine20; R = 0±200")

