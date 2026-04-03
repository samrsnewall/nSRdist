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
%dA = load("../Results/dataT_All1_RLGtrue_DS0p05_Jun2_fitJun24_depthweight.mat");
dA = load("../Results/dataT_All1_RLGtrue_DS0p05_Jun2_fitJul16_depthweight.mat");

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
MLNacc      = NaN(length(dStrus),1);
LNacc      = NaN(length(dStrus),1);
InvGamAcc   = NaN(length(dStrus),1);
GamAcc   = NaN(length(dStrus),1);
BMacc    = NaN(length(dStrus),1);

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
end
dA.summT = table(MeanAgePairsT, MeanSedLength, MeanSedTime, nSR_median,...
    nSR_min, nSR_max, nSR_95lo, nSR_95hi, LNacc, MLNacc, GamAcc, InvGamAcc,...
    BMacc, 'RowNames', dStrusStrings);

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

%% All Bchron and RSR500 samplings, histogram with 68th percentile bars
nsubs = 4;
numruns = length(dA.d.S1.BSampIR.OneRunDatas);
 figure;
% subplot(nsubs,1,1)
% histogram(log(BM.hist),'BinEdges', logBinEdges, 'FaceColor', 'b', 'FaceAlpha', 0.1)
% hold on
% plot(BM.lnSR.x, 1500*BM.lnSR.px, 'LineWidth', 2)
% box on
% xlim([-2.5 2.5])
% ylim([0 3500])
% ylabel("cm")
% title("BIGMACS - Lee et al., 2023")
subplot(nsubs,1,2-1)
histogram(log(dA.d.S1.BMode.weightedC),'BinEdges', logBinEdges, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts');
hold on
plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
box on
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BMode")
subplot(nsubs,1,3-1)
A = sort(dA.d.S1.BSampIR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BSamp [500-4000yr]")
subplot(nsubs,1,4-1)
A = sort(dA.d.S1.New0IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
%xlabel("log(nSR)")
title("RSR0")
subplot(nsubs,1,5-1)
A = sort(dA.d.S1.New500IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
xlabel("log(nSR)")
title("RSR500")

%% Plot the many different fits

figure;
hold on
plot(dA.d.S1.New500IR.MLN.lnSR.x, dA.d.S1.New500IR.MLN.lnSR.px, 'k')
plot(dA.d.S1.New500IR.invGam.lnSR.x, dA.d.S1.New500IR.invGam.lnSR.px, 'r')


%% Load data of just A7 core
d1C = load("../Results/dataT_RLGtrue_R200M20_Apr23_A7_fitApr23_depthweight.mat");

figure;
subplot(4,1,1)
histogram(log(d1C.d.S1.BMode.weightedC),'BinEdges', logBinEdges, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts');
box on
xlim([-2.5 2.5])
%ylim([0 100])
ylabel("cm")
title("BMode")
subplot(4,1,2)
A = sort(d1C.d.S1.BSampIR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', mean(A), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
%ylim([0 100])
ylabel("cm")
title("BSamp [500-4000yr]")
subplot(4,1,3)
A = sort(d1C.d.S1.New0IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', mean(A), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
%ylim([0 100])
ylabel("cm")
title("RSR0")
subplot(4,1,4)
A = sort(d1C.d.S1.New500IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', mean(A), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
%ylim([0 100])
ylabel("cm")
xlabel("log(nSR)")
title("RSR500")

%% All Bchron and RSR500 samplings, histogram with 68th percentile bars
nsubs = 3;
numruns = length(dA.d.S1.BSampIR.OneRunDatas);
 figure;
% subplot(nsubs,1,1)
% histogram(log(BM.hist),'BinEdges', logBinEdges, 'FaceColor', 'b', 'FaceAlpha', 0.1)
% hold on
% plot(BM.lnSR.x, 1500*BM.lnSR.px, 'LineWidth', 2)
% box on
% xlim([-2.5 2.5])
% ylim([0 3500])
% ylabel("cm")
% title("BIGMACS - Lee et al., 2023")
subplot(nsubs,1,2-1)
histogram(log(dA.d.S1.BMode.weightedC),'BinEdges', logBinEdges, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts');
hold on
plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
box on
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BMode")
subplot(nsubs,1,3-1)
A = sort(dA.d.S1.BSampIR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BSamp [500-4000yr]")
subplot(nsubs,1,4-1)
A = sort(dA.d.S1.New500IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
xlabel("log(nSR)")
title("RSR500")

%% Load data of just A7 core
d1C = load("../Results/dataT_RLGtrue_R200M20_Apr23_A7_fitApr23_depthweight.mat");

figure;
subplot(3,1,1)
histogram(log(d1C.d.S1.BMode.weightedC),'BinEdges', logBinEdges, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts');
box on
xlim([-2.5 2.5])
%ylim([0 100])
ylabel("cm")
title("BMode")
subplot(3,1,2)
A = sort(d1C.d.S1.BSampIR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', mean(A), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
%ylim([0 100])
ylabel("cm")
title("BSamp [500-4000yr]")

subplot(3,1,3)
A = sort(d1C.d.S1.New500IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', mean(A), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
%ylim([0 100])
ylabel("cm")
xlabel("log(nSR)")
title("RSR500")

%% Plot A7 age depth model probabilistic

% A7Bchron.input = readtable("../BchronFolders/Bchron_RLGtrue_R200M20_May14_DS0p05/Outputs/A7/inputData.txt");
% A7Bchron.theta = readmatrix("../BchronFolders/Bchron_RLGtrue_R200M20_May14_DS0p05/Outputs/A7/theta.csv");


%figure;
%plot(d1C.d.dataT.bchronProb)