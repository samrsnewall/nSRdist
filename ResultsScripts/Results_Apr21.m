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
dA = load("../Results/dataT_RLGtrue_R200M20_Mar4_fitMar31.mat");

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
plotAgeModes(dA.d.S1.chooseLog, dA.d.dataT.ageModes, dA.d.dataT.cores)

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
    MLNacc(i,1) = sum(dStru.MLN.chiStats.h == 0)./numRunsEva;
    InvGamAcc(i,1) = sum(dStru.invGam.chiStats.h == 0)./numRunsEva;
    GamAcc(i,1) = sum(dStru.Gam.chiStats.h == 0)./numRunsEva;
    LNacc(i,1) = sum(dStru.LN.chiStats.h == 0)./numRunsEva;
end
dA.summT = table(MeanAgePairsT, MeanSedLength, MeanSedTime, nSR_median,...
    nSR_min, nSR_max, nSR_95lo, nSR_95hi,  MLNacc, InvGamAcc,...
    'RowNames', dStrusStrings)

%% Test BM vs results of BMode and new data
figure;
subplot(2,1,1)
yyaxis left
histogram(BM.hist)
yyaxis right
plot(BM.nSR.x, BM.nSR.px)
xlim([0 6])
subplot(2,1,2)
yyaxis left
histogram(dA.d.S1.BMode.weightedC)
yyaxis right
plot(dA.d.S1.BMode.MLN.nSR.x, dA.d.S1.BMode.MLN.nSR.px)
xlim([0 6])

fitS.enforceBinSizeLimits =0;
fitS.dispChi2 = false;
chi2_dataVSpdfVEC(log(dA.d.S1.BMode.weightedC), dA.d.S1.BMode.numCpairs, 20, BMpdf, fitS)

%% Test BM vs results of BSamp and new data

%% Test BM vs results of RSR500 and new data


%% RSRsamplings Histogram with median and internal 68th percentile bars for each bin
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

%% All Bchron and RSR500 samplings, histogram with 68th percentile bars
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
A = sort(dA.d.S1.BSampIR.lnSRHistCounts, 1);
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
OneNewMLNpdf.pdfName = "2 Component Mixed Log Norm";
OneNewMLNpdf.numParams = 6;
[OneNewMLNpdf.x, OneNewMLNpdf.px] = px_to_pfx(dA.d.x,OneNewMLN, f_log);
OneNewInvGam = dA.d.S1.New500IR.invGam.nSR.px(:,r);
OneNewInvGampdf.pdfName = "Inverse Gamma";
OneNewInvGampdf.numParams = 2;
[OneNewInvGampdf.x, OneNewInvGampdf.px] = px_to_pfx(dA.d.x,OneNewInvGam, f_log);
fitS.dispChi2 = true;
fitS.chi2MinCountNum = 5;

pdfs = {OneNewMLNpdf; OneNewInvGampdf};

[h, p, chiStat] = chi2_dataVStwopdfVECs(log(dA.d.S1.New500IR.weightedC{r}), dA.d.S1.New500IR.numCpairs(r), 20,pdfs, fitS);
gcf
xlabel("log(nSR)")

%% Compare Newall sampled data and fits to Bchron sampled data and fits

%Find which distribution fits better for each run
chooseMLN_B1 = sum(dA.d.S1.BSampIR.MLN.chiStats.p > dA.d.S1.BSampIR.invGam.chiStats.p)./dA.d.S1.fitS.OneRun.numruns;
chooseInvGam_B1 = 1 - chooseMLN_B1;
acceptMLN_B1 = sum(dA.d.S1.BSampIR.MLN.chiStats.p > 0.05)/dA.d.S1.fitS.OneRun.numruns;
acceptInvGam_B1 = sum(dA.d.S1.BSampIR.invGam.chiStats.p > 0.05)/dA.d.S1.fitS.OneRun.numruns;
acceptNchooseMLN_B1 = sum(dA.d.S1.BSampIR.MLN.chiStats.p > 0.05 & dA.d.S1.BSampIR.MLN.chiStats.p > dA.d.S1.BSampIR.invGam.chiStats.p)/dA.d.S1.fitS.OneRun.numruns;
acceptNchooseInvGam_B1 = sum(dA.d.S1.BSampIR.invGam.chiStats.p > 0.05 & dA.d.S1.BSampIR.MLN.chiStats.p < dA.d.S1.BSampIR.invGam.chiStats.p)/dA.d.S1.fitS.OneRun.numruns;

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
A = sort(dA.d.S1.BSampIR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
yyaxis right
hold on
plot(dA.d.S1.BSampIR.MLN.lnSR.x, dA.d.S1.BSampIR.MLN.lnSR.px,'LineStyle', '-', 'Color', [0,0,1,0.1],'Marker','none' ,'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [0,0,1], 'Marker','none', 'DisplayName', "MLNs: " + num2str(acceptNchooseMLN_B1*100, 2) + "%")
plot(dA.d.S1.BSampIR.invGam.lnSR.x, dA.d.S1.BSampIR.invGam.lnSR.px,'LineStyle', '-', 'Color', [1,0,0,0.1],'Marker','none', 'HandleVisibility', 'off')
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
chooseMLN_B1 = sum(dA.d.S1.BSampIR.MLN.chiStats.p > dA.d.S1.BSampIR.invGam.chiStats.p)./dA.d.S1.fitS.OneRun.numruns;
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
A = sort(dA.d.S1.BSampIR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
yyaxis right
hold on
plot(dA.d.S1.BSampIR.MLN.lnSR.x, dA.d.S1.BSampIR.MLN.lnSR.px,'LineStyle', '-', 'Color', [0,0,1,0.1],'Marker','none' ,'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [0,0,1], 'Marker','none', 'DisplayName', "MLNs: " + num2str(chooseMLN_B1))
plot(dA.d.S1.BSampIR.invGam.lnSR.x, dA.d.S1.BSampIR.invGam.lnSR.px,'LineStyle', '-', 'Color', [1,0,0,0.1],'Marker','none', 'HandleVisibility', 'off')
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
BT = histogram(log(cell2mat(dADT.d.S1.BSampIR.weightedC)), 'BinEdges', logBinEdges, 'FaceColor', [0, 0, 0.8]);
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
BF = histogram(log(cell2mat(dADF.d.S1.BSampIR.weightedC)), 'BinEdges', logBinEdges, 'FaceColor', [0, 0, 0.8]);
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

%% Evaluate unweighted results

dNW = load("../Results/dataT_RLGtrue_R200M20_Mar4_fitApr8_noweight.mat");

%Create summary table of useful information
dNW.d.S1.BMode.weightedC = dNW.d.S1.BMode.weightedC;
dNW.d.S1.BMode.MLN.chiStats = dNW.d.S1.BMode.MLN.chiStats;
dNW.d.S1.BMode.invGam.chiStats = dNW.d.S1.BMode.invGam.chiStats;
dStrus = {dNW.d.S1.BMode, dNW.d.S1.BSampIR, dNW.d.S1.New0IR, dNW.d.S1.New500IR, dNW.d.S1.New1000IR, dNW.d.S1.New1500IR};
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
dNW.summT = table(MeanAgePairsT, MeanSedLength, MeanSedTime, nSR_median,...
    nSR_min, nSR_max, nSR_95lo, nSR_95hi,  MLNacc, InvGamAcc,...
    'RowNames', dStrusStrings)

%% Evaluate unweighted results

dAW = load("../Results/dataT_RLGtrue_R200M20_Mar4_fitApr8_ageweight.mat");

%Create summary table of useful information
dAW.d.S1.BMode.weightedC = dAW.d.S1.BMode.weightedC;
dAW.d.S1.BMode.MLN.chiStats = dAW.d.S1.BMode.MLN.chiStats;
dAW.d.S1.BMode.invGam.chiStats = dAW.d.S1.BMode.invGam.chiStats;
dStrus = {dAW.d.S1.BMode, dAW.d.S1.BSampIR, dAW.d.S1.New0IR, dAW.d.S1.New500IR, dAW.d.S1.New1000IR, dAW.d.S1.New1500IR};
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
dAW.summT = table(MeanAgePairsT, MeanSedLength, MeanSedTime, nSR_median,...
    nSR_min, nSR_max, nSR_95lo, nSR_95hi,  MLNacc, InvGamAcc,...
    'RowNames', dStrusStrings)


