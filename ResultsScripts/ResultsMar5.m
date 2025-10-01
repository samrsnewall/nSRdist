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

%% Load data and fits
dA = load("../Results/dataT_RLGtrue_R200M20_Mar4_fitMar31.mat");

%%
%Create table of useful information
%summT = table('VariableNames',{'AgePairs', 'SedLength', 'SedTime', 'nSRmedian', 'nSRmin', 'nSRmax', 'InvGamRej', 'MLNRej'});

dA.d.S1.BMode.weightedC = dA.d.S1.BMode.weightedC;
dA.d.S1.BMode.MLN.chiStats = dA.d.S1.BMode.MLN.chiStats;
dA.d.S1.BMode.invGam.chiStats = dA.d.S1.BMode.invGam.chiStats;
dStrus = {dA.d.S1.BMode, dA.d.S1.BChIR, dA.d.S1.New0IR, dA.d.S1.New500IR, dA.d.S1.New1000IR, dA.d.S1.New1500IR};
dStrusStrings = ["BMode","BSamp", "RSR0", "RSR500", "RSR1000", "RSR1500"];

MeanAgePairsT  = NaN(length(dStrus),1);
MeanSedLength  = NaN(length(dStrus),1);
MeanSedTime    = NaN(length(dStrus),1);
nSR_median     = NaN(length(dStrus),1);
nSR_min     = NaN(length(dStrus),1);
nSR_max     = NaN(length(dStrus),1);
nSR_95lo     = NaN(length(dStrus),1);
nSR_95hi     = NaN(length(dStrus),1);
MLNacc      = NaN(length(dStrus),1);
InvGamAcc   = NaN(length(dStrus),1);

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
    MLNacc(i,1) = sum(dStru.MLN.chiStats.h == 0);
    InvGamAcc(i,1) = sum(dStru.invGam.chiStats.h == 0);

end
summT = table(MeanAgePairsT, MeanSedLength, MeanSedTime, nSR_median,...
    nSR_min, nSR_max, nSR_95lo, nSR_95hi,  MLNacc, InvGamAcc,...
    'RowNames', dStrusStrings)

%% Figure of BM histogram in nSR, log(nSR) and inv(nSR)
figure;
subplot(3,1,1)
xlim([0 7])
ylim([0 2000])
hold on
ylabel("cm")
histogram(BM.hist, 'BinEdges', binEdges)
xlabel("nSR")
subplot(3,1,2)
hold on
histogram(log(BM.hist), 'BinEdges', logBinEdges)
xlim([-2.5 2.5])
ylim([0 2000])
xlabel("log(nSR)")
ylabel("cm")
subplot(3,1,3)
hold on
histogram(1./(BM.hist), 'BinEdges', invBinEdges)
xlim([0 7])
ylim([0 2000])
xlabel("inv(nSR)")
ylabel("cm")
ax = gca;
ax.XDir = 'reverse';

%% RSRsamplings Histogram with median and internal 95th percentile bars for each bin
numruns = length(dA.d.S1.New0IR.lnSRHistCounts);

figure;
subplot(4,1,1)
A = sort(dA.d.S1.New0IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylims = [0 4000];
ylim(ylims)
ylabel("cm")
title("RSR0")
subplot(4,1,2)
A = sort(dA.d.S1.New500IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim(ylims)
ylabel("cm")
title("RSR500")
subplot(4,1,3)
A = sort(dA.d.S1.New1000IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim(ylims)
ylabel("cm")
title("RSR1000")
subplot(4,1,4)
A = sort(dA.d.S1.New1500IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim(ylims)
ylabel("cm")
title("RSR1500")
xlabel("log(nSR)")

%% RSRsamplings Histogram with median and internal 95th percentile bars for each bin
numruns = length(dA.d.S1.New0IR.lnSRHistCounts);

figure;
subplot(4,1,1)
A = sort(dA.d.S1.New0IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 5000])
ylabel("cm")
title("RSR0")
subplot(4,1,2)
A = sort(dA.d.S1.New500IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 5000])
ylabel("cm")
title("RSR500")
subplot(4,1,3)
A = sort(dA.d.S1.New1000IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 5000])
ylabel("cm")
title("RSR1000")
subplot(4,1,4)
A = sort(dA.d.S1.New1500IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 5000])
ylabel("cm")
title("RSR1500")
xlabel("log(nSR)")

%% All Bchron and RSR500 samplings, histogram with 95 percents

figure;
subplot(4,1,1)
histogram(log(BM.hist),'BinEdges', logBinEdges, 'FaceColor', 'b', 'FaceAlpha', 0.1)
box on
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BIGMACS - Lee et al., 2023")
subplot(4,1,2)
histogram(log(dA.d.S1.BMode.weightedC),'BinEdges', logBinEdges, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts');
box on
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BMode")
subplot(4,1,3)
A = sort(dA.d.S1.BChIR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BSamp")
subplot(4,1,4)
A = sort(dA.d.S1.New500IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
xlabel("log(nSR)")
title("RSR500")

%% Chi-squared GOF test on 1 sampling of Newall500 sampling approach - with evenly spaced bars (except for end bars)
rng(2)
r = randi(numruns, 1);

f_log = @(x) log(x);
OneNewSamp = dA.d.S1.New500IR.weightedC{r};
OneNewMLN = dA.d.S1.New500IR.MLN.nSR.px(:,r);
OneNewMLNpdf.numParams = 6;
[OneNewMLNpdf.x, OneNewMLNpdf.px] = px_to_pfx(dA.d.x,OneNewMLN, f_log);
OneNewInvGam = dA.d.S1.New500IR.invGam.nSR.px(:,r);
OneNewInvGampdf.numParams = 2;
[OneNewInvGampdf.x, OneNewInvGampdf.px] = px_to_pfx(dA.d.x,OneNewInvGam, f_log);
fitS.dispChi2 = true;

[h1, p1, chiStat1, h2, p2, chiStat2] = chi2_dataVStwopdfVECs(log(dA.d.S1.New500IR.weightedC{r}), dA.d.S1.New500IR.numCpairs(r), 20,OneNewMLNpdf , OneNewInvGampdf, fitS);
gcf
xlabel("log(nSR)")


%% Compare Newall sampled data and fits to Bchron sampled data and fits

%Find which distribution fits better for each run
chooseMLN_B1 = sum(dA.d.S1.BChIR.MLN.chiStats.p > dA.d.S1.BChIR.invGam.chiStats.p)./dA.d.S1.fitS.OneRun.numruns;
chooseInvGam_B1 = 1 - chooseMLN_B1;
acceptMLN_B1 = sum(dA.d.S1.BChIR.MLN.chiStats.p > 0.05)/dA.d.S1.fitS.OneRun.numruns;
acceptInvGam_B1 = sum(dA.d.S1.BChIR.invGam.chiStats.p > 0.05)/dA.d.S1.fitS.OneRun.numruns;
acceptNchooseMLN_B1 = sum(dA.d.S1.BChIR.MLN.chiStats.p > 0.05 & dA.d.S1.BChIR.MLN.chiStats.p > dA.d.S1.BChIR.invGam.chiStats.p)/dA.d.S1.fitS.OneRun.numruns;
acceptNchooseInvGam_B1 = sum(dA.d.S1.BChIR.invGam.chiStats.p > 0.05 & dA.d.S1.BChIR.MLN.chiStats.p < dA.d.S1.BChIR.invGam.chiStats.p)/dA.d.S1.fitS.OneRun.numruns;

%Find which distribution fits better for each run
chooseMLN_Newall500 = sum(dA.d.S1.New500IR.MLN.chiStats.p > dA.d.S1.New500IR.invGam.chiStats.p)./dA.d.S1.fitS.OneRun.numruns;
chooseInvGam_Newall500 = 1 - chooseMLN_Newall500;
acceptMLN_Newall500 = sum(dA.d.S1.New500IR.MLN.chiStats.p > 0.05)/dA.d.S1.fitS.OneRun.numruns;
acceptInvGam_Newall500 = sum(dA.d.S1.New500IR.invGam.chiStats.p > 0.05)/dA.d.S1.fitS.OneRun.numruns;
acceptNchooseMLN_Newall500 = sum(dA.d.S1.New500IR.MLN.chiStats.p > 0.05 & dA.d.S1.New500IR.MLN.chiStats.p > dA.d.S1.New500IR.invGam.chiStats.p)/dA.d.S1.fitS.OneRun.numruns;
acceptNchooseInvGam_Newall500 = sum(dA.d.S1.New500IR.invGam.chiStats.p > 0.05 & dA.d.S1.New500IR.MLN.chiStats.p < dA.d.S1.New500IR.invGam.chiStats.p)/dA.d.S1.fitS.OneRun.numruns;

figure()
subplot(2,1,1)
xlim([-3 3])
yyaxis left
hold on
ylabel("cm")
A = sort(dA.d.S1.BChIR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
yyaxis right
hold on
plot(dA.d.S1.BChIR.MLN.lnSR.x, dA.d.S1.BChIR.MLN.lnSR.px,'LineStyle', '-', 'Color', [0,0,1,0.1],'Marker','none' ,'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [0,0,1], 'Marker','none', 'DisplayName', "MLNs: " + num2str(acceptNchooseMLN_B1*100, 2) + "%")
plot(dA.d.S1.BChIR.invGam.lnSR.x, dA.d.S1.BChIR.invGam.lnSR.px,'LineStyle', '-', 'Color', [1,0,0,0.1],'Marker','none', 'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [1,0,0], 'Marker','none', 'DisplayName', "InvGammas: " + num2str(acceptNchooseInvGam_B1*100, 2) + "%")
ylim([0 1.25])
xlabel("log(nSR)")
ylabel("PDF")
legend()
title("BSamp")
subplot(2,1,2)
xlim([-3 3])
yyaxis left
hold on
ylabel("cm")
A = sort(dA.d.S1.New500IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
yyaxis right
hold on
plot(dA.d.S1.New500IR.MLN.lnSR.x, dA.d.S1.New500IR.MLN.lnSR.px,'LineStyle', '-', 'Color', [0,0,1,0.1],'Marker','none' ,'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [0,0,1], 'Marker','none', 'DisplayName', "MLNs: " + num2str(acceptNchooseMLN_Newall500*100, 2) + "%")
plot(dA.d.S1.New500IR.invGam.lnSR.x, dA.d.S1.New500IR.invGam.lnSR.px,'LineStyle', '-', 'Color', [1,0,0,0.1],'Marker','none', 'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [1,0,0], 'Marker','none', 'DisplayName', "InvGammas: " + num2str(acceptNchooseInvGam_Newall500*100, 2) + "%")
ylim([0 1.25])
xlabel("log(nSR)")
ylabel("PDF")
title("RSR500")
legend()

%% figure with fits with RSR500 and RSR1500

%Find which distribution fits better for each run
chooseMLN_Newall500 = sum(dA.d.S1.New500IR.MLN.chiStats.p > dA.d.S1.New500IR.invGam.chiStats.p)./dA.d.S1.fitS.OneRun.numruns;
chooseInvGam_Newall500 = 1 - chooseMLN_Newall500;
acceptMLN_Newall500 = sum(dA.d.S1.New500IR.MLN.chiStats.p > 0.05)/dA.d.S1.fitS.OneRun.numruns;
acceptInvGam_Newall500 = sum(dA.d.S1.New500IR.invGam.chiStats.p > 0.05)/dA.d.S1.fitS.OneRun.numruns;
acceptNchooseMLN_Newall500 = sum(dA.d.S1.New500IR.MLN.chiStats.p > 0.05 & dA.d.S1.New500IR.MLN.chiStats.p > dA.d.S1.New500IR.invGam.chiStats.p)/dA.d.S1.fitS.OneRun.numruns;

%Find which distribution fits better for each run
chooseMLN_Newall1500 = sum(dA.d.S1.New1500IR.MLN.chiStats.p > dA.d.S1.New1500IR.invGam.chiStats.p)./dA.d.S1.fitS.OneRun.numruns;
chooseInvGam_Newall1500 = 1 - chooseMLN_Newall1500;
acceptMLN_Newall1500 = sum(dA.d.S1.New1500IR.MLN.chiStats.p > 0.05)/dA.d.S1.fitS.OneRun.numruns;
acceptInvGam_Newall1500 = sum(dA.d.S1.New1500IR.invGam.chiStats.p > 0.05)/dA.d.S1.fitS.OneRun.numruns;
acceptNchooseMLN_Newall1500 = sum(dA.d.S1.New1500IR.MLN.chiStats.p > 0.05 & dA.d.S1.New1500IR.MLN.chiStats.p > dA.d.S1.New1500IR.invGam.chiStats.p)/dA.d.S1.fitS.OneRun.numruns;


figure()
subplot(2,1,1)
xlim([-3 3])
yyaxis left
hold on
ylabel("cm")
A = sort(dA.d.S1.New500IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
yyaxis right
hold on
plot(dA.d.S1.New500IR.MLN.lnSR.x, dA.d.S1.New500IR.MLN.lnSR.px,'LineStyle', '-', 'Color', [0,0,1,0.1],'Marker','none' ,'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [0,0,1], 'Marker','none', 'DisplayName', "MLNs: " + num2str(chooseMLN_Newall500))
plot(dA.d.S1.New500IR.invGam.lnSR.x, dA.d.S1.New500IR.invGam.lnSR.px,'LineStyle', '-', 'Color', [1,0,0,0.1],'Marker','none', 'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [1,0,0], 'Marker','none', 'DisplayName', "InvGammas: " + num2str(chooseInvGam_Newall500))
ylim([0 1.25])
xlabel("log(nSR)")
ylabel("PDF")
legend()
title("RSR500")
subplot(2,1,2)
xlim([-3 3])
yyaxis left
hold on
ylabel("cm")
A = sort(dA.d.S1.New1500IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
yyaxis right
hold on
plot(dA.d.S1.New1500IR.MLN.lnSR.x, dA.d.S1.New1500IR.MLN.lnSR.px,'LineStyle', '-', 'Color', [0,0,1,0.1],'Marker','none' ,'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [0,0,1], 'Marker','none', 'DisplayName', "MLNs: " + num2str(chooseMLN_Newall1500))
plot(dA.d.S1.New1500IR.invGam.lnSR.x, dA.d.S1.New1500IR.invGam.lnSR.px,'LineStyle', '-', 'Color', [1,0,0,0.1],'Marker','none', 'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [1,0,0], 'Marker','none', 'DisplayName', "InvGammas: " + num2str(chooseInvGam_Newall1500))
ylim([0 1.25])
xlabel("log(nSR)")
ylabel("PDF")
title("RSR1500")
legend()

%% Compare Newall sampled data and fits to Bchron sampled data and fits

%Find which distribution fits better for each run
chooseMLN_B1 = sum(dA.d.S1.BChIR.MLN.chiStats.p > dA.d.S1.BChIR.invGam.chiStats.p)./dA.d.S1.fitS.OneRun.numruns;
chooseInvGam_B1 = 1 - chooseMLN_B1;

%Find which distribution fits better for each run
chooseMLN_Newall500 = sum(dA.d.S1.New500IR.MLN.chiStats.p > dA.d.S1.New500IR.invGam.chiStats.p)./dA.d.S1.fitS.OneRun.numruns;
chooseInvGam_Newall500 = 1 - chooseMLN_Newall500;

figure()
subplot(2,1,1)
xlim([-3 3])
ylim([0 3000])
yyaxis left
hold on
ylabel("cm")
A = sort(dA.d.S1.BChIR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
yyaxis right
hold on
plot(dA.d.S1.BChIR.MLN.lnSR.x, dA.d.S1.BChIR.MLN.lnSR.px,'LineStyle', '-', 'Color', [0,0,1,0.1],'Marker','none' ,'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [0,0,1], 'Marker','none', 'DisplayName', "MLNs: " + num2str(chooseMLN_B1))
plot(dA.d.S1.BChIR.invGam.lnSR.x, dA.d.S1.BChIR.invGam.lnSR.px,'LineStyle', '-', 'Color', [1,0,0,0.1],'Marker','none', 'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [1,0,0], 'Marker','none', 'DisplayName', "InvGammas: " + num2str(chooseInvGam_B1))
ylim([0 1.25])
xlabel("nSR")
ylabel("PDF")
legend()
title("BSamp")
subplot(2,1,2)
xlim([-3 3])
ylim([0 3000])
yyaxis left
hold on
ylabel("cm")
A = sort(dA.d.S1.New500IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
yyaxis right
hold on
plot(dA.d.S1.New500IR.MLN.lnSR.x, dA.d.S1.New500IR.MLN.lnSR.px,'LineStyle', '-', 'Color', [0,0,1,0.1],'Marker','none' ,'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [0,0,1], 'Marker','none', 'DisplayName', "MLNs: " + num2str(chooseMLN_Newall500))
plot(dA.d.S1.New500IR.invGam.lnSR.x, dA.d.S1.New500IR.invGam.lnSR.px,'LineStyle', '-', 'Color', [1,0,0,0.1],'Marker','none', 'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [1,0,0], 'Marker','none', 'DisplayName', "InvGammas: " + num2str(chooseInvGam_Newall500))
ylim([0 1.25])
xlabel("nSR")
ylabel("PDF")
title("RSR500")
legend()

%% Check fits with changing restriction on Newall sampling approach

%Convert all Newall0 pdfs to log space
datasetN0 = dA.d.S1.New0IR.weightedC;
dN0_label = "Newall0";

%Find which distribution fits better for each run
chooseMLN_Newall0 = sum(dA.d.S1.New0IR.MLN.chiStats.p > dA.d.S1.New0IR.invGam.chiStats.p)./dA.d.S1.fitS.OneRun.numruns;
chooseInvGam_Newall0 = 1 - chooseMLN_Newall0;

%Find which distribution fits better for each run
chooseMLN_Newall1000 = sum(dA.d.S1.New1000IR.MLN.chiStats.p > dA.d.S1.New1000IR.invGam.chiStats.p)./dA.d.S1.fitS.OneRun.numruns;
chooseInvGam_Newall1000 = 1 - chooseMLN_Newall1000;

%Find which distribution fits better for each run
chooseMLN_Newall1500 = sum(dA.d.S1.New1500IR.MLN.chiStats.p > dA.d.S1.New1500IR.invGam.chiStats.p)./dA.d.S1.fitS.OneRun.numruns;
chooseInvGam_Newall1500 = 1 - chooseMLN_Newall1500;

% figure()
% subplot(4,1,1)
% xlim([-3 3])
% yyaxis left
% hold on
% ylabel("cm")
% for i = 1:numruns
%     hold on
%     if i ==1 
%     histogram(log(dA.d.S1.New0IR.weightedC{i}), 'BinEdges', logBinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'on', 'DisplayName', 'Histograms')
%     else
%     histogram(log(dA.d.S1.New0IR.weightedC{i}), 'BinEdges', logBinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'off')
%     end   
% end
% yyaxis right
% hold on
% ylim([0 1])
% xlabel("nSR")
% ylabel("PDF")
% % legend()
% title("Newall0 Samplings")
% subplot(4,1,2)
% xlim([-3 3])
% yyaxis left
% hold on
% ylabel("cm")
% for i = 1:numruns
%     hold on
%     if i ==1 
%     histogram(log(dA.d.S1.New500IR.weightedC{i}), 'BinEdges', logBinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'on', 'DisplayName', 'Histograms')
%     else
%     histogram(log(dA.d.S1.New500IR.weightedC{i}), 'BinEdges', logBinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'off')
%     end   
% end
% yyaxis right
% hold on
% ylim([0 1])
% xlabel("nSR")
% ylabel("PDF")
% title("Newall500 Samplings")
% % legend()
% subplot(4,1,3)
% xlim([-3 3])
% yyaxis left
% hold on
% ylabel("cm")
% for i = 1:numruns
%     hold on
%     if i ==1 
%     histogram(log(dA.d.S1.New1000IR.weightedC{i}), 'BinEdges', logBinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'on', 'DisplayName', 'Histograms')
%     else
%     histogram(log(dA.d.S1.New1000IR.weightedC{i}), 'BinEdges', logBinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'off')
%     end   
% end
% yyaxis right
% hold on
% ylim([0 1])
% xlabel("nSR")
% ylabel("PDF")
% title("Newall1000 Samplings")
% % legend()
% subplot(4,1,4)
% xlim([-3 3])
% yyaxis left
% hold on
% ylabel("cm")
% for i = 1:numruns
%     hold on
%     if i ==1 
%     histogram(log(dA.d.S1.New1500IR.weightedC{i}), 'BinEdges', logBinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'on', 'DisplayName', 'Histograms')
%     else
%     histogram(log(dA.d.S1.New1500IR.weightedC{i}), 'BinEdges', logBinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'off')
%     end   
% end
% yyaxis right
% hold on
% ylim([0 1])
% xlabel("nSR")
% ylabel("PDF")
% title("Newall1500 Samplings")
% legend()

%% 
figure()
subplot(4,1,1)
xlim([-3 3])
yyaxis left
hold on
ylabel("cm")
alph = 0.05;
for i = 1:numruns
    hold on
    if i ==1 
    histogram('BinCounts',dA.d.S1.New0IR.lnSRHistCounts(i,:), 'BinEdges', logBinEdges, 'FaceAlpha', alph, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'on', 'DisplayName', 'Histograms')
    else
    histogram('BinCounts',dA.d.S1.New0IR.lnSRHistCounts(i,:), 'BinEdges', logBinEdges, 'FaceAlpha', alph, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'off')
    end   
end
yyaxis right
hold on
ylim([0 1])
ylabel("PDF")
% legend()
title("RSR0")
subplot(4,1,2)
xlim([-3 3])
yyaxis left
hold on
ylabel("cm")
for i = 1:numruns
    hold on
    if i ==1 
    histogram('BinCounts',dA.d.S1.New500IR.lnSRHistCounts(i,:), 'BinEdges', logBinEdges, 'FaceAlpha', alph, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'on', 'DisplayName', 'Histograms')
    else
    histogram('BinCounts',dA.d.S1.New500IR.lnSRHistCounts(i,:), 'BinEdges', logBinEdges, 'FaceAlpha', alph, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'off')
    end   
end
yyaxis right
hold on
ylim([0 1])
ylabel("PDF")
title("RSR500")
% legend()
subplot(4,1,3)
xlim([-3 3])
yyaxis left
hold on
ylabel("cm")
for i = 1:numruns
    hold on
    if i ==1 
    histogram('BinCounts',dA.d.S1.New1000IR.lnSRHistCounts(i,:), 'BinEdges', logBinEdges, 'FaceAlpha', alph, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'on', 'DisplayName', 'Histograms')
    else
    histogram('BinCounts',dA.d.S1.New1000IR.lnSRHistCounts(i,:), 'BinEdges', logBinEdges, 'FaceAlpha', alph, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'off')
    end   
end
yyaxis right
hold on
ylim([0 1])
ylabel("PDF")
title("RSR1000")
% legend()
subplot(4,1,4)
xlim([-3 3])
yyaxis left
hold on
ylabel("cm")
for i = 1:numruns
    hold on
    if i ==1 
    histogram('BinCounts',dA.d.S1.New1500IR.lnSRHistCounts(i,:), 'BinEdges', logBinEdges, 'FaceAlpha', alph, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'on', 'DisplayName', 'Histograms')
    else
    histogram('BinCounts',dA.d.S1.New1500IR.lnSRHistCounts(i,:), 'BinEdges', logBinEdges, 'FaceAlpha', alph, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'off')
    end   
end
yyaxis right
hold on
ylim([0 1])
xlabel("log(nSR)")
ylabel("PDF")
title("RSR1500")

%% Check impact of using single mean SR vs actual mean SR for sampling

% Load data where meanSR is calculated using sampled data
dADT = load("../Results/dataT_PFandLin_R200M20_Feb19_nSRmeanTrue_fitFeb19.mat");

% Load data where meanSR is calculated using mode of Bchron run
dADF = load("../Results/dataT_PFandLin_R200M20_Feb19_nSRmeanFalse_fitFeb19.mat");

%Create binedges in lognSR
logBinEdges = -5:0.1:5;
logBinCenters = (logBinEdges(1:end-1)+logBinEdges(2:end))./2;

figure;
subplot(3,2,1)
NT = histogram(log(cell2mat(dADT.d.S1.New500IR.weightedC)), 'BinEdges', logBinEdges, 'FaceColor', [0.8, 0, 0]);
xlim([-4 4])
ylim([0 1e6])
xlabel("log(nSR)")
ylabel("cm")
title("True, Newall500")
subplot(3,2,2)
BT = histogram(log(cell2mat(dADT.d.S1.BChIR.weightedC)), 'BinEdges', logBinEdges, 'FaceColor', [0, 0, 0.8]);
xlim([-4 4])
ylim([0 1e6])
xlabel("log(nSR)")
ylabel("cm")
title("True, Bchron")
subplot(3,2,3)
NF = histogram(log(cell2mat(dADF.d.S1.New500IR.weightedC)), 'BinEdges', logBinEdges, 'FaceColor', [0.8, 0, 0]);
xlim([-4 4])
ylim([0 1e6])
xlabel("log(nSR)")
ylabel("cm")
title("False, Newall500")
subplot(3,2,4)
BF = histogram(log(cell2mat(dADF.d.S1.BChIR.weightedC)), 'BinEdges', logBinEdges, 'FaceColor', [0, 0, 0.8]);
xlim([-4 4])
ylim([0 1e6])
xlabel("log(nSR)")
ylabel("cm")
title("False, Bchron")
subplot(3,2,5)
plot(logBinCenters, NT.BinCounts - NF.BinCounts, 'r')
ylim([-1e6 1e6])
title("True - False Bincounts")
xlabel("log(nSR)")
subplot(3,2,6)
plot(logBinCenters, BT.BinCounts - BF.BinCounts, 'b')
ylim([-1e6 1e6])
title("True - False Bincounts")
xlabel("log(nSR)")

%% Compare histograms from depth-weighted, age-weighted, unweighted

dNW = load("../Results/dataT_RLGtrue_R200M20_Mar4_fitApr8_noweight.mat");

