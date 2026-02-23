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
[MLN_BIGMACS, ~, BM.gmfit] = fitMixLogNorm(BM.hist, BM.nSR.x, 2, 3, 343);
[BM.lnSR.x, BM.lnSR.px] = px_to_pfx(BM.nSR.x, BM.nSR.px, @log);
BMpdf.numParams = 6;
BMpdf.x = BM.lnSR.x;
BMpdf.px = BM.lnSR.px;

%% Load data and fits
dA = load("../Results/dataT_All1_RLGtrue_DS0p05_Dec9_fit13Feb26_depthweight_400R");

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

%% Plot age modes 
%Find out which cores passed filtering
usedLog = ~isnan(dA.d.dataT.meanSR);

plotAgeModes(dA.d.S1.chooseLog,dA.d.S1.chooseLog, dA.d.dataT.ageModes, dA.d.dataT.cores)

%% %Create table of useful information
dStrus = {dA.d.S1.BMedian, dA.d.S1.BChIR, dA.d.S1.New0IR, dA.d.S1.New500IR, dA.d.S1.New1000IR};
dStrusStrings = ["BMedian", "BSamp", "RSR0", "RSR500", "RSR1000"];

MeanAgePairsT  = NaN(length(dStrus),1);
MeanSedLength  = NaN(length(dStrus),1);
MeanSedTime    = NaN(length(dStrus),1);
nSR_median     = NaN(length(dStrus),1);
nSR_min     = NaN(length(dStrus),1);
nSR_max     = NaN(length(dStrus),1);
nSR_95lo     = NaN(length(dStrus),1);
nSR_95hi     = NaN(length(dStrus),1);

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
    
    MeanAgePairsT(i,1) = mean(dStru.numCpairs);
    allWCsorted = sort(allWC);
    numWC = length(allWC);
    
    nSR_median(i,1) = median(allWC);
    nSR_min(i,1) = min(allWC);
    nSR_max(i,1) = max(allWC);
    nSR_95lo(i,1) = allWCsorted(ceil(0.025*numWC));
    nSR_95hi(i,1) = allWCsorted(ceil(0.975*numWC));
end
dA.summT = table(MeanAgePairsT, MeanSedLength, MeanSedTime, nSR_median,...
    nSR_min, nSR_max, nSR_95lo, nSR_95hi, 'RowNames', dStrusStrings);

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

%% Plot of BIGMACS and BMedian histograms with BIGMACS distribution superimposed
nsubs = 2;
figure
subplot(nsubs,1,2-1)
histogram(log(BM.hist),'BinEdges', logBinEdges, 'FaceColor', 'b', 'FaceAlpha', 0.1)
hold on
plot(BM.lnSR.x, 1500*BM.lnSR.px, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--')
box on
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BIGMACS: var = " + num2str(var(log(BM.hist)), 3))
subplot(nsubs,1,3-1)
histogram(log(dA.d.S1.BMedian.weightedC),'BinEdges', logBinEdges, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts');
hold on
plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--')
box on
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BMedian: var = " + num2str(var(log(dA.d.S1.BMedian.weightedC)), 3))

%% Plot of above with normalized histogram
nsubs = 2;
fA = figure;
subplot(nsubs,1,2-1)
BM.logHC = histcounts(log(BM.hist), logBinEdges);
normF = 1./sum((BM.logHC.*uniquetol(diff(logBinEdges), 1e-5)));
hold on
histogram('BinCounts', BM.logHC, 'BinEdges', logBinEdges, 'FaceColor', 'b', 'FaceAlpha', 0.1)
plot(BM.lnSR.x, BM.lnSR.px.*(1/normF), 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--')
box on
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BIGMACS")
subplot(nsubs,1,3-1)
BMedian_logHC = histcounts(log(dA.d.S1.BMedian.weightedC), logBinEdges)./dA.d.S1.fitS.DeterministicRun.weightInflator;
normF = 1./sum((BMedian_logHC.*uniquetol(diff(logBinEdges), 1e-5)));
histogram('BinCounts', BMedian_logHC,'BinEdges', logBinEdges, 'FaceColor', 'r', 'FaceAlpha', 0.1);
hold on
plot(BM.lnSR.x, BM.lnSR.px*(1/normF), 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--')
box on
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BMedian")
xlabel("log(nSR)")

%saveas(fA, "../Figures/BMvsBMedian_normalizedPDF.jpg")

%% Plot of BIGMACS and BMedian with normalized pdfs and histogram of dt

load("../Results/Lin2014_dts.mat")

%Get weighted histogram of dts for BIGMACS
BM_dt_weightedBC = makeWeightedBinCounts(asr(:,2), asr(:,3), 0:0.1:10);

nsubs = 4;
fA = figure;
subplot(2,2,2-1)
BM.logHC = histcounts(log(BM.hist), logBinEdges);
normF = 1./sum((BM.logHC.*uniquetol(diff(logBinEdges), 1e-5)));
hold on
histogram('BinCounts', BM.logHC, 'BinEdges', logBinEdges, 'FaceColor', 'b', 'FaceAlpha', 0.1)
plot(BM.lnSR.x, BM.lnSR.px.*(1/normF), 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--')
box on
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BIGMACS")
subplot(2,2,2)
histogram('BinCounts', BM_dt_weightedBC, 'BinEdges', 0:0.1:10, 'FaceColor', 'b', 'FaceAlpha', 0.1)
ylabel("cm")
xlabel("dt")
xlim([0 5])
ylim([0 2000])
subplot(2,2,3)
BMedian_logHC = histcounts(log(dA.d.S1.BMedian.weightedC), logBinEdges);
normF = 1./sum((BMedian_logHC.*uniquetol(diff(logBinEdges), 1e-5)));
histogram('BinCounts', BMedian_logHC,'BinEdges', logBinEdges, 'FaceColor', 'r', 'FaceAlpha', 0.1);
hold on
plot(BM.lnSR.x, BM.lnSR.px*(1/normF), 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--')
box on
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BMedian")
xlabel("log(nSR)")
subplot(2, 2, 4)

histogram('BinCounts', dA.d.S1.BMedian.agediffsWC, 'BinEdges', 0:0.1:10, 'FaceColor', 'r', 'FaceAlpha', 0.1)
ylabel("cm")
xlabel("dt")
xlim([0 5])
ylim([0 2000])

%% Plot of pooled samples + fits

AR = dA.d.S1.New500AR;
IR = dA.d.S1.New500IR;

figure;
subplot(2,1,1)
hold on
AR_logHC = histcounts(log(AR.weightedC), logBinEdges);
normF = 1./sum((AR_logHC.*uniquetol(diff(logBinEdges), 1e-5)));
histogram('BinCounts', AR_logHC, 'BinEdges', logBinEdges, 'HandleVisibility', 'off', 'FaceColor', '[0.8 0.8 0.8]')
plot(AR.MLN.lnSR.x, AR.MLN.lnSR.px*(1/normF), 'r', 'LineWidth', 1, 'DisplayName','MLN')
plot(AR.LN.lnSR.x, AR.LN.lnSR.px*(1/normF), 'LineWidth', 1, 'DisplayName','LN')
plot(AR.invGam.lnSR.x, AR.invGam.lnSR.px*(1/normF),  'LineWidth', 1, 'DisplayName','invGam', 'Color', 'b')
plot(AR.Gam.lnSR.x, AR.Gam.lnSR.px*(1/normF),  'LineWidth', 1, 'DisplayName','Gam' )
xlim([-2.5 2.5])
xlabel("log(nSR)")
ylabel("Pooled Counts")
legend()
ylim([0 1.2e5])
ylims = ylim;

%And the individual runs as well
subplot(2,1,2)
hold on;
AR_logHC = histcounts(log(AR.weightedC), logBinEdges);
normF = 1./sum((AR_logHC.*uniquetol(diff(logBinEdges), 1e-5)));
histogram('BinCounts', AR_logHC, 'BinEdges', logBinEdges, 'HandleVisibility', 'off', 'FaceColor', '[0.8 0.8 0.8]')
plot(IR.MLN.lnSR.x, IR.MLN.lnSR.px(:,1).*(1/normF), 'Color', [0 0 0 0.1], 'DisplayName', 'Individual Run MLN')
plot(IR.MLN.lnSR.x, IR.MLN.lnSR.px(:,2:end).*(1/normF), 'Color', [0 0 0 0.1], 'HandleVisibility', 'off')
plot(AR.MLN.lnSR.x, AR.MLN.lnSR.px.*(1/normF), 'Color', "r", 'LineWidth', 1, 'DisplayName', 'MLN')
xlim([-2.5 2.5])
ylim(ylims)
xlabel("log(nSR)")
ylabel("Pooled Counts")
legend()

%% Double check the pooled weighted histograms match the pooled weighted histograms of individual runs

IR_WCpool = [];
for i = 1:length(IR.weightedC)
    IR_WCpool = [IR_WCpool, IR.weightedC{i}];
end

figure;
subplot(2,1,1)
histogram('BinCounts', AR_logHC, 'BinEdges', logBinEdges, 'HandleVisibility', 'off', 'FaceColor', '[0.8 0.8 0.8]')
subplot(2,1,2)
IR_logHC = histcounts(log(IR_WCpool), logBinEdges);
histogram('BinCounts', IR_logHC, 'BinEdges', logBinEdges, 'HandleVisibility', 'off', 'FaceColor', '[0.8 0.8 0.8]')


%% All Bchron and RSR500 samplings, histogram with 68th percentile bars
nsubs = 3;
numruns = length(dA.d.S1.BChIR.OneRunDatas);
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
histogram(log(dA.d.S1.BMedian.weightedC),'BinEdges', logBinEdges, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts');
hold on
plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
box on
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BMedian [500-4000yr]")
subplot(nsubs,1,3-1)
BSamp_hists = sort(dA.d.S1.BChIR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', BSamp_hists(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BSamp [500-4000yr]")
subplot(nsubs,1,4-1)
MCSamp_hists = sort(dA.d.S1.New500IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', MCSamp_hists(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
xlabel("log(nSR)")
title("MCSamp")


%% All sampling approaches histogram of median with 68th percentile bars + dt histogram
nsubs = 4;
numruns = length(dA.d.S1.BChIR.OneRunDatas);
figure;
subplot(4,2,1)
A = sort(dA.d.S1.BChIR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
%plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
xlabel("log(nSR)")
ylabel("cm")
title("BSamp [500-4000yr]")
subplot(4,2,2)
B = sort(dA.d.S1.BChIR.agediffsWbC', 1);
hold on
box on
histogram('BinCounts', B(numruns*0.5, :), 'BinEdges', 0:0.1:10, 'FaceAlpha', 0.1)
errorbar(0.05:0.1:9.95, B(numruns*0.5, :),(B(numruns*0.5, :)-B(numruns*0.16, :)),(B(numruns*0.84, :)-B(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
xlim([0 6])
ylim([0 2000])
xlabel("dt (kyr)")
ylabel("cm")
subplot(4,2,3)
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
xlabel("log(nSR)")
ylabel("cm")
%xlabel("log(nSR)")
title("RSR0")
subplot(4,2,4)
B = sort(dA.d.S1.New0IR.agediffsWbC', 1);
hold on
box on
histogram('BinCounts', B(numruns*0.5, :), 'BinEdges', 0:0.1:10, 'FaceAlpha', 0.1)
errorbar(0.05:0.1:9.95, B(numruns*0.5, :),(B(numruns*0.5, :)-B(numruns*0.16, :)),(B(numruns*0.84, :)-B(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
xlim([0 6])
ylim([0 2000])
xlabel("dt (kyr)")
ylabel("cm")


subplot(4,2,5)
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
xlabel("log(nSR)")
ylabel("cm")
title("RSR500")
subplot(4,2,6)
B = sort(dA.d.S1.New500IR.agediffsWbC', 1);
hold on
box on
histogram('BinCounts', B(numruns*0.5, :), 'BinEdges', 0:0.1:10, 'FaceAlpha', 0.1)
errorbar(0.05:0.1:9.95, B(numruns*0.5, :),(B(numruns*0.5, :)-B(numruns*0.16, :)),(B(numruns*0.84, :)-B(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
xlim([0 6])
ylim([0 2000])
xlabel("dt (kyr)")
ylabel("cm")


subplot(4,2,7)
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
xlabel("log(nSR)")
ylabel("cm")
title("RSR1000")
subplot(4,2,8)
B = sort(dA.d.S1.New1000IR.agediffsWbC', 1);
hold on
box on
histogram('BinCounts', B(numruns*0.5, :), 'BinEdges', 0:0.1:10, 'FaceAlpha', 0.1)
errorbar(0.05:0.1:9.95, B(numruns*0.5, :),(B(numruns*0.5, :)-B(numruns*0.16, :)),(B(numruns*0.84, :)-B(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
xlim([0 6])
ylim([0 2000])
xlabel("dt (kyr)")
ylabel("cm")

%% All sampling approaches histogram pooled
nsubs = 4;
numruns = length(dA.d.S1.BChIR.OneRunDatas);
figure;
subplot(4,2,1)
A = sort(dA.d.S1.BChIR.lnSRHistCounts, 1);
box on
histogram('BinCounts', sum(A,1), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
xlim([-2.5 2.5])
%ylim([0 3500])
xlabel("log(nSR)")
ylabel("cm")
title("BSamp [500-4000yr]")
subplot(4,2,2)
B = sort(dA.d.S1.BChIR.agediffsWbC', 1);
hold on
box on
histogram('BinCounts', sum(B,1), 'BinEdges', 0:0.1:10, 'FaceAlpha', 0.1)
xlim([0 6])
%ylim([0 2000])
xlabel("dt (kyr)")
ylabel("cm")
subplot(4,2,3)
A = sort(dA.d.S1.New0IR.lnSRHistCounts, 1);
box on
histogram('BinCounts', sum(A,1), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
xlim([-2.5 2.5])
%ylim([0 3500])
xlabel("log(nSR)")
ylabel("cm")
%xlabel("log(nSR)")
title("RSR0")
subplot(4,2,4)
B = sort(dA.d.S1.New0IR.agediffsWbC', 1);
hold on
box on
histogram('BinCounts', sum(B,1), 'BinEdges', 0:0.1:10, 'FaceAlpha', 0.1)
xlim([0 6])
%ylim([0 2000])
xlabel("dt (kyr)")
ylabel("cm")


subplot(4,2,5)
A = sort(dA.d.S1.New500IR.lnSRHistCounts, 1);
box on
histogram('BinCounts', sum(A,1), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
xlim([-2.5 2.5])
%ylim([0 3500])
xlabel("log(nSR)")
ylabel("cm")
title("RSR500")
subplot(4,2,6)
B = sort(dA.d.S1.New500IR.agediffsWbC', 1);
hold on
box on
histogram('BinCounts', sum(B,1), 'BinEdges', 0:0.1:10, 'FaceAlpha', 0.1)
xlim([0 6])
%ylim([0 2000])
xlabel("dt (kyr)")
ylabel("cm")

subplot(4,2,7)
A = sort(dA.d.S1.New1000IR.lnSRHistCounts, 1);
box on
histogram('BinCounts', sum(A,1), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
xlim([-2.5 2.5])
%ylim([0 3500])
xlabel("log(nSR)")
ylabel("cm")
title("RSR1000")
subplot(4,2,8)
B = sort(dA.d.S1.New1000IR.agediffsWbC', 1);
hold on
box on
histogram('BinCounts', sum(B,1), 'BinEdges', 0:0.1:10, 'FaceAlpha', 0.1)
xlim([0 6])
%ylim([0 2000])
xlabel("dt (kyr)")
ylabel("cm")

%% Compare the variance of each approach
commonxlim = [0 5];

figure;
subplot(4,1,1)
xline(var(BM.hist))
xlim(commonxlim)
subplot(4,1,2)
xline(var(dA.d.S1.BMedian.weightedC))
xlim(commonxlim)
subplot(4,1,3)
histogram(cellfun(@var, dA.d.S1.BChIR.weightedC), 0:0.2:5)
xlim(commonxlim)
subplot(4,1,4)
histogram(cellfun(@var, dA.d.S1.New500IR.weightedC),  0:0.2:5)
xlim(commonxlim)

%% BIC (taeheefix) of each method
%This is the BIC (with the correction to the likelihood) for each 

BIC_BMedian = [dA.d.S1.BMedian.LN.fits{1}.BICtaeheefix, dA.d.S1.BMedian.MLN.fits{1}.BICtaeheefix, dA.d.S1.BMedian.Gam.fits{1}.BICtaeheefix, dA.d.S1.BMedian.invGam.fits{1}.BICtaeheefix]
BIC_BChAR = [dA.d.S1.BChAR.LN.fitInfo.BICtaeheefix, dA.d.S1.BChAR.MLN.fitInfo.BICtaeheefix, dA.d.S1.BChAR.Gam.fitInfo.BICtaeheefix, dA.d.S1.BChAR.invGam.fitInfo.BICtaeheefix]
BIC_New0AR = [dA.d.S1.New0AR.LN.fitInfo.BICtaeheefix, dA.d.S1.New0AR.MLN.fitInfo.BICtaeheefix, dA.d.S1.New0AR.Gam.fitInfo.BICtaeheefix, dA.d.S1.New0AR.invGam.fitInfo.BICtaeheefix]
BIC_New500AR = [dA.d.S1.New500AR.LN.fitInfo.BICtaeheefix, dA.d.S1.New500AR.MLN.fitInfo.BICtaeheefix, dA.d.S1.New500AR.Gam.fitInfo.BICtaeheefix, dA.d.S1.New500AR.invGam.fitInfo.BICtaeheefix]
BIC_New1000AR = [dA.d.S1.New1000AR.LN.fitInfo.BICtaeheefix, dA.d.S1.New1000AR.MLN.fitInfo.BICtaeheefix, dA.d.S1.New1000AR.Gam.fitInfo.BICtaeheefix, dA.d.S1.New1000AR.invGam.fitInfo.BICtaeheefix]

%% Inverse Gamma Parameters for each method

invGam_BMedian = [dA.d.S1.BMedian.invGam.fits{1}.alpha, dA.d.S1.BMedian.invGam.fits{1}.beta];
invGam_BChAR = [dA.d.S1.BChAR.invGam.fitInfo.alpha, dA.d.S1.BChAR.invGam.fitInfo.beta];
invGam_New0AR = [dA.d.S1.New0AR.invGam.fitInfo.alpha, dA.d.S1.New0AR.invGam.fitInfo.beta];
invGam_New500AR = [dA.d.S1.New500AR.invGam.fitInfo.alpha, dA.d.S1.New500AR.invGam.fitInfo.beta];
invGam_New1000AR =[dA.d.S1.New1000AR.invGam.fitInfo.alpha, dA.d.S1.New1000AR.invGam.fitInfo.beta];

invGam_parameters = [invGam_BMedian; invGam_BChAR; invGam_New0AR; invGam_New500AR; invGam_New1000AR]

%% Gamma parameters for each method

Gam_BMedian = [dA.d.S1.BMedian.Gam.fits{1}.alpha, dA.d.S1.BMedian.Gam.fits{1}.beta];
Gam_BChAR = [dA.d.S1.BChAR.Gam.fitInfo.alpha, dA.d.S1.BChAR.Gam.fitInfo.beta];
Gam_New0AR = [dA.d.S1.New0AR.Gam.fitInfo.alpha, dA.d.S1.New0AR.Gam.fitInfo.beta];
Gam_New500AR = [dA.d.S1.New500AR.Gam.fitInfo.alpha, dA.d.S1.New500AR.Gam.fitInfo.beta];
Gam_New1000AR =[dA.d.S1.New1000AR.Gam.fitInfo.alpha, dA.d.S1.New1000AR.Gam.fitInfo.beta];

Gam_parameters = [Gam_BMedian; Gam_BChAR; Gam_New0AR; Gam_New500AR; Gam_New1000AR]

%% LN parameters for each method

LN_BMedian = [dA.d.S1.BMedian.LN.fits{1}.mu, dA.d.S1.BMedian.LN.fits{1}.Sigma];
LN_BChAR = [dA.d.S1.BChAR.LN.fitInfo.mu, dA.d.S1.BChAR.LN.fitInfo.Sigma];
LN_New0AR = [dA.d.S1.New0AR.LN.fitInfo.mu, dA.d.S1.New0AR.LN.fitInfo.Sigma];
LN_New500AR = [dA.d.S1.New500AR.LN.fitInfo.mu, dA.d.S1.New500AR.LN.fitInfo.Sigma];
LN_New1000AR =[dA.d.S1.New1000AR.LN.fitInfo.mu, dA.d.S1.New1000AR.LN.fitInfo.Sigma];

LN_parameters = [LN_BMedian; LN_BChAR; LN_New0AR; LN_New500AR; LN_New1000AR]