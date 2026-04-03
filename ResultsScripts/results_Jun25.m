%% Add necessary paths
addpath("../Functions/")

%% Set up histogram settings
%Create binedges
maxbinedge   = 20;
bw           = 0.1; %bin width
binEdges     =  0:bw:maxbinedge; %nSR
invBinEdges  =  0:bw:maxbinedge; %invnSR
logBinEdges  = -5:bw:5;          %lognSR
logBinCenters = (logBinEdges(1:end-1)+logBinEdges(2:end)).*0.5;

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
BMpdf.numParams = 6;
BMpdf.x = BM.lnSR.x;
BMpdf.px = BM.lnSR.px;

%% Load data and fits
dA = load("../Results/dataT_All1_RLGtrue_DS0p05_Jun2_fitAug22_depthweight_1000R.mat");
%dA = load("../Results/dataT_All1_RLGtrue_DS0p05_Jun2_fitJul4_noweight.mat");
%% Take data of interest
%----- Find corenames, lats, longs, depths of cores included
dataTbl = dA.d.dataT(dA.d.S1.chooseLog, :);
%----- Make map with locations denoted as red squares
figure;
subplot(2,2,[1 2])
worldmap('World')
setm(gca, 'mapprojection', 'robinson')
geoshow('landareas.shp', 'FaceColor','[0.7 0.7 0.7]', 'EdgeColor', '[0.7 0.7 0.7]')
load coastlines coastlat coastlon
plotm(coastlat, coastlon, 'Color', 'k')
hold on
plotm(dataTbl.lats, dataTbl.longs, 'rs')
%Plot histograms of Mean SR, Depths, Resolution by Age, Resolution by Depth

subplot(2,2,3)
histogram(dataTbl.depths./1000, 0:0.25:6, 'FaceColor', 'k')
xlim([0 6])
xlabel("Depth (km)")
ylabel("Counts")
xticks(0:0.5:6)
title("Core Depths")

subplot(2,2,4)
histogram(dataTbl.meanSR, 0:2:90, 'FaceColor','k')
xlabel('Mean SR (cm/kyr)')
ylabel('Counts')
xlim([0 90])
title("Cores' Mean SR")

%% Plot age modes (Need to remind myself how age modes were calculated, are they through bchron? the large gaps between ages are kept in)
%Find out which cores passed filtering
usedLog = ~isnan(dA.d.dataT.meanSR);

plotAgeModes(dA.d.S1.chooseLog,dA.d.S1.chooseLog, dA.d.dataT.ageModes, dA.d.dataT.cores)

%% %Create table of useful information
dA.d.S1.BMode.weightedC = dA.d.S1.BMode.weightedC;
dA.d.S1.BMode.MLN.chiStats = dA.d.S1.BMode.MLN.chiStats;
dA.d.S1.BMode.invGam.chiStats = dA.d.S1.BMode.invGam.chiStats;
dStrus = {dA.d.S1.BMode, dA.d.S1.BSampIR, dA.d.S1.New0IR, dA.d.S1.New500IR, dA.d.S1.New1000IR, dA.d.S1.New1500IR};
dStrusStrings = ["BMode","BSamp", "RSR0", "RSR500", "RSR1000", "RSR1500"];
numRunsEva = length(dA.d.S1.BSampIR.MLN.chiStats.h);

MeanAgePairsT  = NaN(length(dStrus),1);
MeanSedLength  = NaN(length(dStrus),1);
MeanSedTime    = NaN(length(dStrus),1);
nSR_median     = NaN(length(dStrus),1);
nSR_min     = NaN(length(dStrus),1);
nSR_max     = NaN(length(dStrus),1);
nSR_95lo     = NaN(length(dStrus),1);
nSR_95hi     = NaN(length(dStrus),1);
BMacc    = NaN(length(dStrus),1);
MLNacc      = NaN(length(dStrus),1);
LNacc      = NaN(length(dStrus),1);
InvGamAcc   = NaN(length(dStrus),1);
GamAcc   = NaN(length(dStrus),1);
MLNbest      = NaN(length(dStrus),1);
LNbest      = NaN(length(dStrus),1);
InvGambest   = NaN(length(dStrus),1);
Gambest   = NaN(length(dStrus),1);
MLNaccbest      = NaN(length(dStrus),1);
LNaccbest      = NaN(length(dStrus),1);
InvGamaccbest   = NaN(length(dStrus),1);
Gamaccbest   = NaN(length(dStrus),1);

for i = 1:length(dStrus)
    dStru = dStrus{i};
    if isa(dStru.weightedC, 'cell')
        allWC = cell2mat(dStru.weightedC);
            MeanSedLength(i,1) = mean(dStru.sedLength);
    MeanSedTime(i,1) = mean(dStru.sedTimeSpan);
    else
        allWC = dStru.weightedC;
        MeanSedLength(i,1) = mean(dStru.sedLength);
        MeanSedTime(i,1)   = mean(dStru.sedTimeSpan);
    end
    numRunsEva = length(dStru.MLN.chiStats.h);
    
    MeanAgePairsT(i,1) = mean(dStru.numCpairs);
    allWCsorted = sort(allWC);
    numWC = length(allWC);
    
    nSR_median(i,1) = median(allWC);
    nSR_min(i,1) = min(allWC);
    nSR_max(i,1) = max(allWC);
    nSR_95lo(i,1) = allWCsorted(ceil(0.025*numWC));
    nSR_95hi(i,1) = allWCsorted(ceil(0.975*numWC));
    MLNacc(i,1) = sum(dStru.MLN.chiStats.h == 0)./numRunsEva;
    InvGamAcc(i,1) = sum(dStru.invGam.chiStats.h == 0)./numRunsEva;
    GamAcc(i,1) = sum(dStru.Gam.chiStats.h == 0)./numRunsEva;
    LNacc(i,1) = sum(dStru.LN.chiStats.h == 0)./numRunsEva;
    BMacc(i,1) = sum(dStru.chiStatTvsBM.h == 0)./numRunsEva; 
    
    MLNp = dStru.MLN.chiStats.p;
    LNp = dStru.LN.chiStats.p;
    Gamp = dStru.Gam.chiStats.p;
    InvGamp = dStru.invGam.chiStats.p;
    MLNbest(i,1) = sum(MLNp > max([LNp, Gamp, InvGamp], [], 2))./numRunsEva;
    LNbest(i,1) = sum(LNp > max([MLNp, Gamp, InvGamp], [], 2))./numRunsEva;
    InvGambest(i,1) = sum(InvGamp > max([MLNp, Gamp, LNp], [], 2))./numRunsEva;
    Gambest(i,1) = sum(Gamp > max([MLNp, InvGamp, LNp], [], 2))./numRunsEva;

    MLNaccbest(i,1) = sum(dStru.MLN.chiStats.h == 0 & MLNp > max([LNp, Gamp, InvGamp], [], 2))./numRunsEva;
    LNaccbest(i,1) = sum(dStru.LN.chiStats.h == 0 & LNp > max([MLNp, Gamp, InvGamp], [], 2))./numRunsEva;
    InvGamaccbest(i,1) = sum(dStru.invGam.chiStats.h == 0 & InvGamp > max([MLNp, Gamp, LNp], [], 2))./numRunsEva;
    Gamaccbest(i,1) = sum(dStru.Gam.chiStats.h == 0 & Gamp > max([MLNp, InvGamp, LNp], [], 2))./numRunsEva;
end
dA.summT = table(MeanAgePairsT, MeanSedLength, MeanSedTime, nSR_median,...
    nSR_min, nSR_max, nSR_95lo, nSR_95hi, BMacc, LNacc, MLNacc, GamAcc, InvGamAcc,...
    MLNbest, LNbest, InvGambest, Gambest, MLNaccbest,...
    LNaccbest, InvGamaccbest, Gamaccbest,'RowNames', dStrusStrings);

%% BIGMACS Histogram and fit
figure;
histogram(log(BM.hist),'BinEdges', logBinEdges, 'FaceColor', 'b', 'FaceAlpha', 0.1)
hold on
plot(BM.lnSR.x, 1500*BM.lnSR.px, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--')
box on
xlim([-2.5 2.5])
ylim([0 2000])
ylabel("cm")
xlabel("log(nSR)")
title("BIGMACS - Lee et al., 2023")

%% Perform chi2GOF on BIGMACS histogram and distribution
dA.d.S1.fitS.dispChi2 = true;
[BM.chi2.h, BM.chi2.p, BM.chi2.chistat] = chi2gof_vsMLN(BM.gmfit, log(BM.hist), 343, dA.d.S1.fitS);

%% Recreate chi2GOF figure for BMode data vs BM distribution
dA.d.S1.fitS.dispChi2 = true;
[h, p, chiStat1Run_BM] = chi2gof_vsMLN(BM.gmfit, log(dA.d.S1.BMode.weightedC), dA.d.S1.BMode.numCpairs, dA.d.S1.fitS);

%% Recreate chi2GOF figure for BMode data vs 4 fits

%Make a column cell-vector that holds all pdfs to test
pdfs = {dA.d.S1.BMode.MLN.lnSR; dA.d.S1.BMode.LN.lnSR; dA.d.S1.BMode.Gam.lnSR; dA.d.S1.BMode.invGam.lnSR};
%pdfs = {d.S1.BMode.MLN.lnSR; d.S1.BMode.invGam.lnSR};
h = NaN(1,1);
p = NaN(1,1);
[h, p, chiStat] = chi2_dataVStwopdfVECs(log(dA.d.S1.BMode.weightedC), dA.d.S1.BMode.numCpairs, 20, pdfs, dA.d.S1.fitS);

%% BIGMACS and BMode histogram of median with 68th percentile bars
nsubs = 2;
numruns = length(dA.d.S1.BSampIR.OneRunDatas);
 figure;
subplot(nsubs,1,1)
histogram(log(BM.hist),'BinEdges', logBinEdges, 'FaceColor', 'b', 'FaceAlpha', 0.1, "Normalization", 'pdf', "HandleVisibility", "off")
hold on
plot(BM.lnSR.x, BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', "DisplayName","BIGMACS Distribution")
box on
xlim([-2.5 2.5])
%ylim([0 3500])
ylabel("PDF")
title("BIGMACS - Lee et al., 2023")
subplot(nsubs,1,2)
histogram(log(dA.d.S1.BMode.weightedC),'BinEdges', logBinEdges,"Normalization", 'pdf', 'FaceColor', 'r', 'FaceAlpha', 0.1, "HandleVisibility", "off");
hold on
plot(BM.lnSR.x, BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', "DisplayName","BIGMACS Distribution")
box on
xlim([-2.5 2.5])
%ylim([0 3500])
ylabel("PDF")
title("BMode")
xlabel("log(nSR)")
legend()

%% All sampling approaches histogram of median with 68th percentile bars
nsubs = 5;
numruns = length(dA.d.S1.BSampIR.OneRunDatas);
figure;
subplot(nsubs,1,1)
A = sort(dA.d.S1.BSampIR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
%plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BSamp [500-4000yr]")
subplot(nsubs,1,2)
A = sort(dA.d.S1.New0IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
%plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
%xlabel("log(nSR)")
title("RSR0")
subplot(nsubs,1,3)
A = sort(dA.d.S1.New500IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
%plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("RSR500")
subplot(nsubs,1,4)
A = sort(dA.d.S1.New1000IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
%plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("RSR1000")
subplot(nsubs,1,5)
A = sort(dA.d.S1.New1500IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
%plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("RSR1500")
xlabel("log(nSR)")

%% All sampling approaches, histogram of mean with 68th percentile bars
nsubs = 5;
numruns = length(dA.d.S1.BSampIR.OneRunDatas);
figure;
subplot(nsubs,1,1)
A = sort(dA.d.S1.BSampIR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', mean(A,1), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
%plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BSamp [500-4000yr]")
subplot(nsubs,1,2)
A = sort(dA.d.S1.New0IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', mean(A,1), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
%plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
%xlabel("log(nSR)")
title("RSR0")
subplot(nsubs,1,3)
A = sort(dA.d.S1.New500IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', mean(A,1), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
%plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("RSR500")
subplot(nsubs,1,4)
A = sort(dA.d.S1.New1000IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', mean(A,1), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
%plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("RSR1000")
subplot(nsubs,1,5)
A = sort(dA.d.S1.New1500IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', mean(A,1), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
%plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("RSR1500")
xlabel("log(nSR)")

%% Plot the many different fits of BMode Approach

AtoPlot = dA.d.S1.BMode;

figure;
subplot(2,2,1)
hold on
histogram(log(dA.d.S1.BMode.weightedC),'BinEdges', logBinEdges,"Normalization", 'pdf', 'FaceColor', 'r', 'FaceAlpha', 0.1, "HandleVisibility", "off");
plot(AtoPlot.MLN.lnSR.x, AtoPlot.MLN.lnSR.px, 'Color', 'k')
xlim([-2.5 2.5])
ylim([0 1.25])
xlabel("log(nSR)")
ylabel("pdf")
title("Mixed Lognormal")
box("on")
subplot(2,2,2)
hold on
histogram(log(dA.d.S1.BMode.weightedC),'BinEdges', logBinEdges,"Normalization", 'pdf', 'FaceColor', 'r', 'FaceAlpha', 0.1, "HandleVisibility", "off");
plot(AtoPlot.invGam.lnSR.x, AtoPlot.invGam.lnSR.px, 'Color', 'k')
xlim([-2.5 2.5])
ylim([0 1.25])
xlabel("log(nSR)")
ylabel("pdf")
title("Inverse Gamma")
box("on")
subplot(2,2,3)
hold on
histogram(log(dA.d.S1.BMode.weightedC),'BinEdges', logBinEdges,"Normalization", 'pdf', 'FaceColor', 'r', 'FaceAlpha', 0.1, "HandleVisibility", "off");
plot(AtoPlot.LN.lnSR.x, AtoPlot.LN.lnSR.px, 'Color', 'k')
xlim([-2.5 2.5])
xlabel("log(nSR)")
ylabel("pdf")
title("Lognormal")
box("on")
ylim([0 1.25])
subplot(2,2,4)
hold on
histogram(log(dA.d.S1.BMode.weightedC),'BinEdges', logBinEdges,"Normalization", 'pdf', 'FaceColor', 'r', 'FaceAlpha', 0.1, "HandleVisibility", "off");
plot(AtoPlot.Gam.lnSR.x, AtoPlot.Gam.lnSR.px, 'Color', 'k')
xlim([-2.5 2.5])
ylim([0 1.25])
xlabel("log(nSR)")
ylabel("pdf")
title("Gamma")
box("on")


%% Plot the many different fits of an approach

AtoPlot = dA.d.S1.BMode;

figure;
subplot(2,2,1)
hold on
plot(AtoPlot.MLN.lnSR.x, AtoPlot.MLN.lnSR.px, 'Color', [0 0 0 0.05])
plot(BM.lnSR.x, BM.lnSR.px, 'r','LineWidth', 1)
xlim([-2.5 2.5])
ylim([0 1])
xlabel("log(nSR)")
ylabel("pdf")
title("Mixed Lognormal")
box("on")
subplot(2,2,2)
plot(AtoPlot.invGam.lnSR.x, AtoPlot.invGam.lnSR.px, 'Color', [0 0 0 0.05])
xlim([-2.5 2.5])
ylim([0 1])
xlabel("log(nSR)")
ylabel("pdf")
title("Inverse Gamma")
box("on")
subplot(2,2,3)
plot(AtoPlot.LN.lnSR.x, AtoPlot.LN.lnSR.px, 'Color', [0 0 0 0.05])
xlim([-2.5 2.5])
xlabel("log(nSR)")
ylabel("pdf")
title("Lognormal")
box("on")
ylim([0 1])
subplot(2,2,4)
plot(AtoPlot.Gam.lnSR.x, AtoPlot.Gam.lnSR.px, 'Color', [0 0 0 0.05])
xlim([-2.5 2.5])
ylim([0 1])
xlabel("log(nSR)")
ylabel("pdf")
title("Gamma")
box("on")

%% Plot the many different fits of an approach vs fit to all samples

AtoPlotIR = dA.d.S1.New1000IR;
AtoPlotAR = dA.d.S1.New1000AR;

figure;
subplot(2,2,1)
hold on
plot(AtoPlotIR.MLN.lnSR.x, AtoPlotIR.MLN.lnSR.px, 'Color', [0 0 0 0.05])
plot(BM.lnSR.x, BM.lnSR.px, 'r','LineWidth', 1)
xlim([-2.5 2.5])
ylim([0 1])
xlabel("log(nSR)")
ylabel("pdf")
title("Mixed Lognormal")
box("on")
subplot(2,2,2)
plot(AtoPlotIR.invGam.lnSR.x, AtoPlotIR.invGam.lnSR.px, 'Color', [0 0 0 0.05])
xlim([-2.5 2.5])
ylim([0 1])
xlabel("log(nSR)")
ylabel("pdf")
title("Inverse Gamma")
box("on")
subplot(2,2,3)
plot(AtoPlotIR.LN.lnSR.x, AtoPlotIR.LN.lnSR.px, 'Color', [0 0 0 0.05])
xlim([-2.5 2.5])
xlabel("log(nSR)")
ylabel("pdf")
title("Lognormal")
box("on")
ylim([0 1])
subplot(2,2,4)
plot(AtoPlotIR.Gam.lnSR.x, AtoPlotIR.Gam.lnSR.px, 'Color', [0 0 0 0.05])
xlim([-2.5 2.5])
ylim([0 1])
xlabel("log(nSR)")
ylabel("pdf")
title("Gamma")
box("on")


subplot(2,2,1)
hold on
plot(AtoPlotAR.MLN.lnSR.x, AtoPlotAR.MLN.lnSR.px, 'Color', 'g', 'LineWidth', 1)
xlim([-2.5 2.5])
ylim([0 1])
xlabel("log(nSR)")
ylabel("pdf")
title("Mixed Lognormal")
box("on")
subplot(2,2,2)
hold on
plot(AtoPlotAR.invGam.lnSR.x, AtoPlotAR.invGam.lnSR.px, 'Color', 'g', 'LineWidth', 1)
xlim([-2.5 2.5])
ylim([0 1])
xlabel("log(nSR)")
ylabel("pdf")
title("Inverse Gamma")
box("on")
subplot(2,2,3)
hold on
plot(AtoPlotAR.LN.lnSR.x, AtoPlotAR.LN.lnSR.px, 'Color', 'g', 'LineWidth', 1)
xlim([-2.5 2.5])
xlabel("log(nSR)")
ylabel("pdf")
title("Lognormal")
box("on")
ylim([0 1])
subplot(2,2,4)
hold on
plot(AtoPlotAR.Gam.lnSR.x, AtoPlotAR.Gam.lnSR.px, 'Color', 'g', 'LineWidth', 1)
xlim([-2.5 2.5])
ylim([0 1])
xlabel("log(nSR)")
ylabel("pdf")
title("Gamma")
box("on")

figure;
subplot(2,2,1)
hold on
plot(AtoPlotIR.MLN.lnSR.mu, AtoPlotIR.MLN.lnSR.var, 'kx', 'LineStyle', "none")
plot(AtoPlotAR.MLN.lnSR.mu, AtoPlotAR.MLN.lnSR.var, 'go', "LineStyle", "none", "MarkerFaceColor", "g")
xlabel("mean")
ylabel("var")
title("MLN")
subplot(2,2,3)
hold on
plot(AtoPlotIR.LN.lnSR.mu, AtoPlotIR.LN.lnSR.var, 'kx', 'LineStyle', "none")
plot(AtoPlotAR.LN.lnSR.mu, AtoPlotAR.LN.lnSR.var, 'go', "LineStyle", "none", "MarkerFaceColor", "g")
xlabel("mean")
ylabel("var")
title("Lognormal")
subplot(2,2,2)
hold on
plot(AtoPlotIR.invGam.lnSR.mu, AtoPlotIR.invGam.lnSR.var, 'kx', 'LineStyle', "none")
plot(AtoPlotAR.invGam.lnSR.mu, AtoPlotAR.invGam.lnSR.var, 'go', "LineStyle", "none", "MarkerFaceColor", "g")
xlabel("mean")
ylabel("var")
title("Inverse Gamma")
subplot(2,2,4)
hold on
plot(AtoPlotIR.Gam.lnSR.mu, AtoPlotIR.Gam.lnSR.var, 'kx', 'LineStyle', "none")
plot(AtoPlotAR.Gam.lnSR.mu, AtoPlotAR.Gam.lnSR.var, 'go', "LineStyle", "none", "MarkerFaceColor", "g")
xlabel("mean")
ylabel("var")
title("Gamma")

%% Plot histogram of all combined results

nsubs = 5;
numruns = length(dA.d.S1.BSampIR.OneRunDatas);
figure;
subplot(nsubs,1,1)
yyaxis("left")
box on
histogram(log(cell2mat(dA.d.S1.BSampIR.weightedC)), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
ylabel("cm")
yyaxis("right")
plot(BM.lnSR.x, BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', "DisplayName", "BIGMACS")
xlim([-2.5 2.5])
%ylim([0 3500])
ylabel("pdf")
title("BSamp")
legend()
subplot(nsubs,1,2)

hold on
box on
histogram(log(cell2mat(dA.d.S1.New0IR.weightedC)), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
%plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
xlim([-2.5 2.5])
%ylim([0 3500])
ylabel("cm")
%xlabel("log(nSR)")
title("RSR0")
subplot(nsubs,1,3)

hold on
box on
histogram(log(cell2mat(dA.d.S1.New500IR.weightedC)), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
%plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
xlim([-2.5 2.5])
%ylim([0 3500])
ylabel("cm")
title("RSR500")
subplot(nsubs,1,4)

hold on
box on
histogram(log(cell2mat(dA.d.S1.New1000IR.weightedC)), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
%plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
xlim([-2.5 2.5])
%ylim([0 3500])
ylabel("cm")
title("RSR1000")
subplot(nsubs,1,5)
hold on
box on
histogram(log(cell2mat(dA.d.S1.New1500IR.weightedC)), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
%plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
xlim([-2.5 2.5])
%ylim([0 3500])
ylabel("cm")
title("RSR1500")
xlabel("log(nSR)")



%%
figure;
subplot(2,2,1)
hold on
plot(AtoPlotIR.MLN.nSR.x, AtoPlotIR.MLN.nSR.px, 'Color', [0 0 0 0.05])
plot(BM.nSR.x, BM.nSR.px, 'r','LineWidth', 1)
xlim([0 6])
ylim([0 1])
xlabel("log(nSR)")
ylabel("pdf")
title("Mixed Lognormal")
box("on")
subplot(2,2,2)
plot(AtoPlotIR.invGam.nSR.x, AtoPlotIR.invGam.nSR.px, 'Color', [0 0 0 0.05])
xlim([0 6])
ylim([0 1])
xlabel("log(nSR)")
ylabel("pdf")
title("Inverse Gamma")
box("on")
subplot(2,2,3)
plot(AtoPlotIR.LN.nSR.x, AtoPlotIR.LN.nSR.px, 'Color', [0 0 0 0.05])
xlim([0 6])
xlabel("log(nSR)")
ylabel("pdf")
title("Lognormal")
box("on")
ylim([0 1])
subplot(2,2,4)
plot(AtoPlotIR.Gam.nSR.x, AtoPlotIR.Gam.nSR.px, 'Color', [0 0 0 0.05])
xlim([0 6])
ylim([0 1])
xlabel("log(nSR)")
ylabel("pdf")
title("Gamma")
box("on")


subplot(2,2,1)
hold on
plot(AtoPlotAR.MLN.nSR.x, AtoPlotAR.MLN.nSR.px, 'Color', 'g', 'LineWidth', 1)
xlim([0 6])
ylim([0 1])
xlabel("log(nSR)")
ylabel("pdf")
title("Mixed Lognormal")
box("on")
subplot(2,2,2)
hold on
plot(AtoPlotAR.invGam.nSR.x, AtoPlotAR.invGam.nSR.px, 'Color', 'g', 'LineWidth', 1)
xlim([0 6])
ylim([0 1])
xlabel("log(nSR)")
ylabel("pdf")
title("Inverse Gamma")
box("on")
subplot(2,2,3)
hold on
plot(AtoPlotAR.LN.nSR.x, AtoPlotAR.LN.nSR.px, 'Color', 'g', 'LineWidth', 1)
xlim([0 6])
xlabel("log(nSR)")
ylabel("pdf")
title("Lognormal")
box("on")
ylim([0 1])
subplot(2,2,4)
hold on
plot(AtoPlotAR.Gam.nSR.x, AtoPlotAR.Gam.nSR.px, 'Color', 'g', 'LineWidth', 1)
xlim([0 6])
ylim([0 1])
xlabel("log(nSR)")
ylabel("pdf")
title("Gamma")
box("on")
