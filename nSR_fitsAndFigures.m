%% Load in desired results file
%load("Results/dataT_highRes_Sep12th24.mat")
load("Results/dataT_planktonicF50_Nov1.mat");
addpath('Functions')


%% ------ Define Subsets of interest
% high SR and low SR cores (separated by 8cm/kyr following Lin et al., 2014)
depth1000Log = dataT.depths > 1000;
lowSRCoresLog   = dataT.meanSR<= 8 & depth1000Log;
highSRCoresLog  = dataT.meanSR >8 & depth1000Log;
allCoresLog     = ~isnan(dataT.meanSR) & depth1000Log;

%% Metadata figures
metadataLog = allCoresLog;
% outputMetadataAndSummaryFigures(metadataLog,dataT)

%% Plot comparison figures 
metadataLog1 = highSRCoresLog;
Log1Colour = "r";
metadataLog2 = lowSRCoresLog;
Log2Colour = "b";
% outputMetadataAndSummaryFiguresComparison(metadataLog1, metadataLog2, dataT, dataT, Log1Colour, Log2Colour, "SR > 8cm/kyr", "SR < 8cm/kyr")

%% Set up useful values
% Get x values of interest
lognorm_BIGMACS = readtable("lognormal_BIGMACS.txt");                       % (use x values currently used in BIGMACS)
x = lognorm_BIGMACS.Var1';
numruns2sample = 400;
%hist_BIGMACS = load("../BIGMACSdata/Lin2014_sedrateratio_cm_wo_NaN.txt");

%% MIX LOG NORMAL FIGURES

%% Compare all runs vs individual runs result for a single set up
%numruns2sample = 400;
desiredLog = highSRCoresLog;
desiredRestriction = dataT.nSRcounts1000;
respectiveString = "HighSR1000";


 [mixLogAllruns] = plotSRandResHistograms(desiredRestriction, x, desiredLog, true, 3, 1, 2, 0, respectiveString,true);
[mixLogNorm1Run, c95up, c95down] = SingleRunLogNorms(desiredRestriction, desiredLog, numruns2sample, x,2, true, 3, 4, 0.001);

figure;
hold on
plot(x, mixLogNorm1Run, 'Color', [0.2 0.2 0.2 0.1], 'LineStyle', '-')
plot(x, [c95up; c95down], '--r', 'LineWidth', 2)
plot(x, mixLogAllruns(:,2),  '-r', 'LineWidth',2)
plot(x, lognorm_BIGMACS.Var2, '--k', 'LineWidth',2)
title(respectiveString)
xlim([0 6])
ylim([0 1.5])
xlabel("nSR")

%% 

[highSRAllruns, highSRhistcounts] = plotSRandResHistograms(dataT.nSRcounts1000, x, highSRCoresLog, true, 3, 1, 2, 0, "All Cores",true);
[lowSRAllruns, lowSRhistcounts] = plotSRandResHistograms(dataT.nSRcounts1000, x, lowSRCoresLog, true, 3, 10, 2, 0, "All Cores",true);

figure;
hold on;
plot(x, highSRAllruns(:,2),  '-r', 'LineWidth',2, "DisplayName", "high")
plot(x, lowSRAllruns(:,2),  '-b', 'LineWidth',2, "DisplayName", "low")
plot(x, lognorm_BIGMACS.Var2, 'k--','LineWidth',2, "DisplayName", "BIGMACS")
%title("Atlantic vs Non-Atlantic")
xlim([0 6])
ylim([0 1.5])
xlabel("nSR")
legend()

%%

figure;
subplot(3,1,1)
hold on
plot(x, lognorm_BIGMACS.Var2, 'k--','LineWidth',2, "DisplayName", "BIGMACS")
xlim([0 6])
ylim([0 1])
subplot(3,1,2)
hold on
histogram(highSRhistcounts, "BinEdges", [0:0.1:10, 1000], "Normalization", "pdf", "FaceColor",'k', "FaceAlpha",0.2)
plot(x, highSRAllruns(:,2),  '-r', 'LineWidth',2, "DisplayName", "high")
xlim([0 6])
ylim([0 1])
subplot(3,1,3)
hold on
histogram(lowSRhistcounts, "BinEdges", [0:0.1:10, 1000], "Normalization", "pdf", "FaceColor",'k', "FaceAlpha",0.2)
plot(x, lowSRAllruns(:,2),  '-b', 'LineWidth',2, "DisplayName", "low")
%title("Atlantic vs Non-Atlantic")
xlim([0 6])
ylim([0 1])
xlabel("nSR")
%legend()


%% Compare first half (alphabetically) of highSR cores to other half of highSR cores

desiredLog = highSRCoresLog;
desiredRestriction = dataT.nSRcounts1000;
total = sum(highSRCoresLog);
halftotal = ceil(total/2);
firsthalfLog = cumsum(desiredLog) <= halftotal;
secondhalfLog = cumsum(desiredLog) > halftotal;
desiredLog1 = desiredLog & firsthalfLog;
desiredLog2 = desiredLog & secondhalfLog;

[mLN_1halfLog] = plotSRandResHistograms(desiredRestriction, x, desiredLog1, true, 3, 1, 2, 0, "All Cores",true);
[mLN_2halfLog] = plotSRandResHistograms(desiredRestriction, x, desiredLog2, true, 3, 1, 2, 0, "All Cores",true);

figure;
hold on

plot(x, mLN_1halfLog(:,2),  '-r', 'LineWidth',2)
plot(x, mLN_2halfLog(:,2),  '-b', 'LineWidth',2)
plot(x, lognorm_BIGMACS.Var2, 'k--','LineWidth',2, "DisplayName", "BIGMACS")
title("Alphabetical First Half vs Second Half")
xlim([0 6])
ylim([0 1.5])
xlabel("nSR")

outputMetadataAndSummaryFiguresComparison(desiredLog1, desiredLog2, dataT, dataT, Log1Colour, Log2Colour, "First Half Alphabetical", "Second Half Alphabetical")

%First alphabetical half is mostly in Atlantic, second alphabetical half is
%mostly in Indonesian Archipelago

%% Test how Atlantic cores compare to other cores
% 
% for i = 1:length(ocean);
%     if contains(ocean{i}, "Atlantic");
%         AtlanticLog(i,1) = 1; 
%     else
%         AtlanticLog(i,1) = 0;
%     end
% end
% 
% desiredLog1 = desiredLog & AtlanticLog;
% desiredLog2 = desiredLog & ~AtlanticLog;
% 
% [mLN_Atlantic] = plotSRandResHistograms(desiredRestriction, x, desiredLog1, true, 3, 1, 2, 0, "All Cores",true);
% [mLN_notAtlantic] = plotSRandResHistograms(desiredRestriction, x, desiredLog2, true, 3, 1, 2, 0, "All Cores",true);
% 
% figure;
% hold on;
% plot(x, mLN_Atlantic(:,2),  '-r', 'LineWidth',2, "DisplayName", "Atlantic")
% plot(x, mLN_notAtlantic(:,2),  '-b', 'LineWidth',2, "DisplayName", "Other")
% plot(x, lognorm_BIGMACS.Var2, 'k--','LineWidth',2, "DisplayName", "BIGMACS")
% title("Atlantic vs Non-Atlantic")
% xlim([0 6])
% ylim([0 1.5])
% xlabel("nSR")
% legend()
% 
% outputMetadataAndSummaryFiguresComparison(desiredLog1, desiredLog2, dataT, dataT, Log1Colour, Log2Colour, "Atlantic", "Not Atlantic")


%% Test how the deepest half of the cores compare to the shallower half

depth1to2kmLog = dataT.depths > 1000 & dataT.depths < 2000;

desiredLog1 = desiredLog & depth1to2kmLog;
desiredLog2 = desiredLog & ~depth1to2kmLog;

[mLN_shallower2000] = plotSRandResHistograms(desiredRestriction, x,desiredLog1, true, 3, 1, 2, 0, "All Cores",true);
[mLN_deeper2000] = plotSRandResHistograms(desiredRestriction, x, desiredLog2, true, 3, 1, 2, 0, "All Cores",true);

outputMetadataAndSummaryFiguresComparison(desiredLog1, desiredLog2, dataT,dataT, Log1Colour, Log2Colour, "Depth < 2km", "Depth > 2km")



 %% Plot weighted vs non-weighted results
 % [SR_MixLogNorm1Run0_noweight] = SingleRunLogNorms(nSRcounts, allCoresLog, numruns, x, 0);
 % [SR_MixLogNorm1Run0_weight] = SingleRunLogNorms(nSRcounts, allCoresLog, numruns, x, 1);
 % [allCores_mixLog0_noweight] = plotSRandResHistograms(nSRcounts, x, agediffs, num14cpairs, allCoresLog, 0, 2, 0, "All Cores",0);
 % [allCores_mixLog0_weight] = plotSRandResHistograms(nSRcounts, x, agediffs, num14cpairs, allCoresLog, 1, 2, 0, "All Cores",0);

% figure;
% subplot(2,1,1)
% hold on
% plot(x, SR_MixLogNorm1Run0_noweight, 'Color', [0.2 0.2 0.2 0.1], 'LineStyle', '-')
% plot(x, allCores_mixLog0_noweight(:,2), '-r', 'LineWidth',  2)
% xlim([0 6])
% ylim([0 1])
% title("no weighting")
% 
% subplot(2,1,2)
% hold on
% plot(x, SR_MixLogNorm1Run0_weight, 'Color', [0.2 0.2 0.2 0.1], 'LineStyle', '-')
% plot(x, allCores_mixLog0_weight(:,2), '-k', 'LineWidth',  2)
% xlim([0 6])
% ylim([0 1])
% title("no weighting")
% 
% figure;
% subplot(2,1,1)
% hold on
% plot(x, allCores_mixLog0_noweight(:,2), '-r', 'LineWidth',  2)
% plot(x, allCores_mixLog0_weight(:,2), '-k', 'LineWidth',  2)
% xlim([0 6])
% ylim([0 1])
% subplot(2,1,2)
% hold on
% plot(x, SR_MixLogNorm1Run0_noweight, 'Color', [1 0 0 0.1], 'LineStyle', '-')
% plot(x, SR_MixLogNorm1Run0_weight, 'Color', [0.2 0.2 0.2 0.1], 'LineStyle', '-')
% xlim([0 6])
% ylim([0 1])

%% Compare weighting in Linspace with weighting in Logspace
% [mixLog_logrounding] = plotSRandResHistograms_sensitivity(nSRcounts500, x, agediffs, num14cpairs, allCoresLog, 1, 3 ,1, 2, 0, "Log Round 3dp",1);
% [mixLog_linrounding] = plotSRandResHistograms(nSRcounts500, x, agediffs, num14cpairs, allCoresLog, 1, 3, 1, 2, 0, "Linear Round 3dp",1);
% 
% figure;
% hold on
% plot(x, mixLog_logrounding(:,2),  '-r', 'LineWidth',2, 'DisplayName', "LogRounding")
% plot(x, mixLog_linrounding(:,2),  '-k', 'LineWidth',2, "DisplayName", "LinRounding")
% plot(x, lognorm_BIGMACS.Var2, '--k', 'LineWidth',2, "DisplayName", "BIGMACS")
% title("Rounding nSR in linspace vs logspace")
% xlim([0 6])
% ylim([0 1])
% legend()
% xlabel("nSR")

%% Test stability of certain set up
% numtests = 20;
% setupMixLogNorms = NaN(length(x), numtests);
% 
% for i = 1:numtests
%     setupMixLogNormHOLDER = plotSRandResHistograms(nSRcounts500, x, agediffs, num14cpairs, A7Log, true, 3, 1, 2, 0, "Linear Round 3dp",false);
%     setupMixLogNorms(:,i) = setupMixLogNormHOLDER(:,2);
% end
% 
% figure;
% plot(x, setupMixLogNorms,'Color', [0.2 0.2 0.2 0.1] )
% xlabel("nSR")
% ylabel("pdf")
% xlim([0 6])

%% Extra figures 

%% Plot high SR 500 against low SR 500

% figure;
% hold on
% plot(x, SR_MixLogNorm1Run500_high, 'Color', [1 0 0 0.1], 'LineStyle', '-')
% plot(x, SR_MixLogNorm1Run500_low, 'Color', [0 0 1 0.1], 'LineStyle', '-')
% xlim([0 6])
% ylim([0 1.5])
% xlabel("nSR")
% title("High SR vs Low SR 500y")
% 

%% Plot all 1 run plots, with combined mixed Log Norm and Bigmacs Mixed Log Norm on top
numruns = 400;
regVal = 0;

agediffsBinEdges = 0:500:10000;

% [allCores_mixLog0, allCores0_histData, allCores0_diffData] = plotSRandResHistograms(dataT.nSRcounts, x,  allCoresLog, true, 3, 1000, 2, regVal, "All Cores",true);
% [allCores_mixLog500, allCores500_histData, allCores500_diffData] = plotSRandResHistograms(dataT.nSRcounts500, x, allCoresLog, true, 3, 1000, 2, regVal, "All Cores",true);
[allCores_mixLog1000, allCores1000_histData, allCores1000_diffData] = plotSRandResHistograms(dataT.nSRcounts1000, x,  allCoresLog,  true, 3, 1000, 2, regVal, "All Cores",true);
[allCores_mixLog1500, allCores1500_histData, allCores1500_diffData] = plotSRandResHistograms(dataT.nSRcounts1500, x,  allCoresLog,  true, 3, 1000, 2, regVal, "All Cores",true);
[allCores_mixLog2000, allCores2000_histData, allCores2000_diffData] = plotSRandResHistograms(dataT.nSRcounts2000, x,allCoresLog,  true, 3, 1000, 2, regVal, "All Cores",true);

% [highSRCores_mixLog0, highSRCores0_histData, highSRCores0_diffData] = plotSRandResHistograms(dataT.nSRcounts, x, highSRCoresLog,  true, 3, 1, 2, regVal, "All Cores",true);
% [highSRCores_mixLog500, highSRCores500_histData, highSRCores500_diffData] = plotSRandResHistograms(dataT.nSRcounts500, x, highSRCoresLog,  true, 3, 1, 2, regVal, "All Cores",true);
[highSRCores_mixLog1000, highSRCores1000_histData, highSRCores1000_diffData] = plotSRandResHistograms(dataT.nSRcounts1000, x, highSRCoresLog,  true, 3, 1, 2, regVal, "All Cores",true);
[highSRCores_mixLog1500, highSRCores1500_histData, highSRCores1500_diffData] = plotSRandResHistograms(dataT.nSRcounts1500, x, highSRCoresLog,  true, 3, 1, 2, regVal, "All Cores",true);
% 
% [lowSRCores_mixLog0, lowSRCores0_histData, lowSRCores0_diffData] = plotSRandResHistograms(dataT.nSRcounts, x, lowSRCoresLog, true, 3, 10, 2, regVal, "All Cores",false);
% [lowSRCores_mixLog500, lowSRCores500_histData, lowSRCores500_diffData] = plotSRandResHistograms(dataT.nSRcounts500, x, lowSRCoresLog,  true, 3, 10, 2, regVal, "All Cores",false);
[lowSRCores_mixLog1000, lowSRCores1000_histData, lowSRCores1000_diffData] = plotSRandResHistograms(dataT.nSRcounts1000, x, lowSRCoresLog,  true, 3, 10, 2, regVal, "All Cores",false);
[lowSRCores_mixLog1500, lowSRCores1500_histData, lowSRCores1500_diffData] = plotSRandResHistograms(dataT.nSRcounts1500, x, lowSRCoresLog,  true, 3, 10, 2, regVal, "All Cores",false);

 % %Compute single run mix log normals for all cores, with 0, 500, 1000 yr
 % %restrictions
 % [SR_MixLogNorm1Run0,   SingleRun0_95pctup,    SingleRun0_95pctdown]    = SingleRunLogNorms(nSRcounts,     allCoresLog, numruns, x, true, 3, 1000, regVal);
 % [SR_MixLogNorm1Run500, SingleRun500_95pctup,  SingleRun500_95pctdown]  = SingleRunLogNorms(nSRcounts500,  allCoresLog, numruns, x, true, 3, 1000, regVal);
  [SR_MixLogNorm1Run1000,SingleRun1000_95pctup, SingleRun1000_95pctdown] = SingleRunLogNorms(dataT.nSRcounts1000, allCoresLog, numruns, x, true, 3, 1000, regVal);
 % 
 % [SR_MixLogNorm1Run0_low,    SingleRun0_low_95pctup,    SingleRun0_low_95pctdown]    = SingleRunLogNorms(nSRcounts,     lowSRCoresLog, numruns, x, true, 3, 1000, regVal);
 % [SR_MixLogNorm1Run500_low,  SingleRun500_low_95pctup,  SingleRun500_low_95pctdown]  = SingleRunLogNorms(nSRcounts500,  lowSRCoresLog, numruns, x, true, 3, 1000, 0.0001);
 [SR_MixLogNorm1Run1000_low, SingleRun1000_low_95pctup, SingleRun1000_low_95pctdown] = SingleRunLogNorms(dataT.nSRcounts1000, lowSRCoresLog, numruns, x, true, 3, 1000, 0.0001);
 % 
 % [SR_MixLogNorm1Run0_high,    SingleRun0_high_95pctup,    SingleRun0_high_95pctdown]    = SingleRunLogNorms(nSRcounts,     highSRCoresLog, numruns, x, true, 3, 1000, regVal);
 % [SR_MixLogNorm1Run500_high,  SingleRun500_high_95pctup,  SingleRun500_high_95pctdown]  = SingleRunLogNorms(nSRcounts500,  highSRCoresLog, numruns, x, true, 3, 1000, regVal);
 [SR_MixLogNorm1Run1000_high, SingleRun1000_high_95pctup, SingleRun1000_high_95pctdown] = SingleRunLogNorms(dataT.nSRcounts1000, highSRCoresLog, numruns, x, true, 3, 1000, regVal);

 ylimits = [0 1.5];

%  %% Plot all 1 run plots, with combined mixed Log Norm and Bigmacs Mixed Log Norm on top
% numruns = 400;
% regVal = 0;
% 
% 
% [fitS.allCoresR0.nSR.oneMLN] = plotSRandResHistograms(nSRcounts, x, agediffs, num14cpairs, allCoresLog, true, 3, 1, 2, regVal, "All Cores",false);
% [fitS.allCoresR500.nSR.oneMLN, fitS.allCoresR500.nSR.weightedCounts] = plotSRandResHistograms(nSRcounts500, x, agediffs, num14cpairs, allCoresLog, true, 3, 1, 2, regVal, "All Cores",false);
% [fitS.allCoresR1000.nSR.oneMLN] = plotSRandResHistograms(nSRcounts1000, x, agediffs, num14cpairs, allCoresLog,  true, 3, 1, 2, regVal, "All Cores",false);
% 
% [fitS.highSRCoresR0.nSR.oneMLN] = plotSRandResHistograms(nSRcounts, x, agediffs, num14cpairs, highSRCoresLog,  true, 3, 1, 2, regVal, "All Cores",false);
% [fitS.highSRCoresR500.nSR.oneMLN, fitS.highSRCoresR500.nSR.weightedCounts] = plotSRandResHistograms(nSRcounts500, x, agediffs, num14cpairs, highSRCoresLog,  true, 3, 1, 2, regVal, "All Cores",false);
% [fitS.highSRCoresR1000.nSR.oneMLN] = plotSRandResHistograms(nSRcounts1000, x, agediffs, num14cpairs, highSRCoresLog,  true, 3, 1, 2, regVal, "All Cores",false);
% 
% [fitS.lowSRCoresR0.nSR.oneMLN] = plotSRandResHistograms(nSRcounts, x, agediffs, num14cpairs, lowSRCoresLog, true, 3, 10, 2, regVal, "All Cores",false);
% [fitS.lowSRCoresR500.nSR.oneMLN, fitS.lowSRCoresR500.nSR.weightedCounts] = plotSRandResHistograms(nSRcounts500, x, agediffs, num14cpairs, lowSRCoresLog,  true, 3, 10, 2, regVal, "All Cores",false);
% [fitS.lowSRCoresR1000.nSR.oneMLN] = plotSRandResHistograms(nSRcounts1000, x, agediffs, num14cpairs, lowSRCoresLog,  true, 3, 10, 2, regVal, "All Cores",false);
% 
%  %Compute single run mix log normals for all cores, with 0, 500, 1000 yr
%  %restrictions
%  [fitS.allCoresR0.nSR.multiMLNs,   SingleRun0_95pctup,    SingleRun0_95pctdown]    = SingleRunLogNorms(nSRcounts,     allCoresLog, numruns, x, true, 3, 1000, regVal);
%  [fitS.allCoresR500.nSR.multiMLNs, SingleRun500_95pctup,  SingleRun500_95pctdown]  = SingleRunLogNorms(nSRcounts500,  allCoresLog, numruns, x, true, 3, 1000, regVal);
%  [fitS.allCoresR1000.nSR.multiMLNs,SingleRun1000_95pctup, SingleRun1000_95pctdown] = SingleRunLogNorms(nSRcounts1000, allCoresLog, numruns, x, true, 3, 1000, regVal);
% 
%  [fitS.lowSRCoresR0.nSR.multiMLNs,    SingleRun0_low_95pctup,    SingleRun0_low_95pctdown]    = SingleRunLogNorms(nSRcounts,     lowSRCoresLog, numruns, x, true, 3, 1000, regVal);
%  [fitS.lowSRCoresR500.nSR.multiMLNs,  SingleRun500_low_95pctup,  SingleRun500_low_95pctdown]  = SingleRunLogNorms(nSRcounts500,  lowSRCoresLog, numruns, x, true, 3, 1000, 0.0001);
%  [fitS.lowSRCoresR1000.nSR.multiMLNs, SingleRun1000_low_95pctup, SingleRun1000_low_95pctdown] = SingleRunLogNorms(nSRcounts1000, lowSRCoresLog, numruns, x, true, 3, 1000, regVal);
% 
%  [fitS.highSRCoresR0.nSR.multiMLNs,    SingleRun0_high_95pctup,    SingleRun0_high_95pctdown]    = SingleRunLogNorms(nSRcounts,     highSRCoresLog, numruns, x, true, 3, 1000, regVal);
%  [fitS.highSRCoresR500.nSR.multiMLNs,  SingleRun500_high_95pctup,  SingleRun500_high_95pctdown]  = SingleRunLogNorms(nSRcounts500,  highSRCoresLog, numruns, x, true, 3, 1000, regVal);
%  [fitS.highSRCoresR1000.nSR.multiMLNs, SingleRun1000_high_95pctup, SingleRun1000_high_95pctdown] = SingleRunLogNorms(nSRcounts1000, highSRCoresLog, numruns, x, true, 3, 1000, regVal);
% 
%  

%% All Cores histograms and res histograms
figure;
subplot(3,2,1)
histogram(allCores0_histData, "BinEdges", [0:0.1:10, 1000], "Normalization", "pdf", "FaceColor",'k', "FaceAlpha",0.2)
xlabel("nSR")
xlim([0 6])
subplot(3,2,2)
histogram("BinCounts", allCores0_diffData, "BinEdges", agediffsBinEdges, "FaceColor",'k', "FaceAlpha",0.2)
xlabel("Age Diff (yr)")
ylim([0 4000])

subplot(3,2,3)
histogram(allCores1000_histData, "BinEdges", [0:0.1:10, 1000], "Normalization", "pdf", "FaceColor",'k', "FaceAlpha",0.2)
xlabel("nSR")
xlim([0 6])
subplot(3,2,4)
histogram("BinCounts", allCores1000_diffData, "BinEdges", agediffsBinEdges, "FaceColor",'k', "FaceAlpha",0.2)
xlabel("Age Diff (yr)")
ylim([0 4000])

subplot(3,2,5)
histogram(allCores2000_histData, "BinEdges", [0:0.1:10, 1000], "Normalization", "pdf", "FaceColor",'k', "FaceAlpha",0.2)
xlabel("nSR")
xlim([0 6])
subplot(3,2,6)
histogram("BinCounts", allCores2000_diffData, "BinEdges", agediffsBinEdges, "FaceColor",'k', "FaceAlpha",0.2)
xlabel("Age Diff (yr)")
ylim([0 4000])

%% lowSR Cores histograms and res histograms
figure;
subplot(3,2,1)
histogram(lowSRCores0_histData, "BinEdges", [0:0.1:10, 1000], "Normalization", "pdf", "FaceColor",'k', "FaceAlpha",0.2)
xlabel("nSR")
xlim([0 6])
subplot(3,2,2)
histogram("BinCounts", lowSRCores0_diffData, "BinEdges", agediffsBinEdges, "FaceColor",'k', "FaceAlpha",0.2)
xlabel("Age Diff (yr)")
ylim([0 500])

subplot(3,2,3)
histogram(lowSRCores500_histData, "BinEdges", [0:0.1:10, 1000], "Normalization", "pdf", "FaceColor",'k', "FaceAlpha",0.2)
xlabel("nSR")
xlim([0 6])
subplot(3,2,4)
histogram("BinCounts",lowSRCores500_diffData, "BinEdges", agediffsBinEdges, "FaceColor",'k', "FaceAlpha",0.2)
xlabel("Age Diff (yr)")
ylim([0 500])

subplot(3,2,5)
histogram(lowSRCores1000_histData, "BinEdges", [0:0.1:10, 1000], "Normalization", "pdf", "FaceColor",'k', "FaceAlpha",0.2)
xlabel("nSR")
xlim([0 6])
subplot(3,2,6)
histogram("BinCounts",lowSRCores1000_diffData, "BinEdges", agediffsBinEdges, "FaceColor",'k', "FaceAlpha",0.2)
xlabel("Age Diff (yr)")
ylim([0 500])

 %% Plot all mix log norm results

 ylimits = [0 1.5];


figure;
subplot(3,3,1)
hold on
plot(x, SR_MixLogNorm1Run0, 'Color', [0.2 0.2 0.2 0.1], 'LineStyle', '-')
plot(x, [SingleRun0_95pctup; SingleRun0_95pctdown], '--r', 'LineWidth', 2)
plot(x, allCores_mixLog0(:,2),  '-k', 'LineWidth',2)
plot(x, lognorm_BIGMACS.Var2, '--k', 'LineWidth',2)
title("0 yr restriction")
xlim([0 6])
ylim(ylimits)
xlabel("nSR")
subplot(3,3,2)
hold on
plot(x, SR_MixLogNorm1Run500, 'Color', [0.2 0.2 0.2 0.1], 'LineStyle', '-')
plot(x, [SingleRun500_95pctup; SingleRun500_95pctdown], '--r', 'LineWidth', 2)
plot(x, allCores_mixLog500(:,2),  '-k', 'LineWidth',2)
plot(x, lognorm_BIGMACS.Var2, '--k', 'LineWidth',2)
title("500 yr restriction")
xlim([0 6])
ylim(ylimits)
xlabel("nSR")
subplot(3,3,3)
hold on
plot(x, SR_MixLogNorm1Run1000, 'Color', [0.2 0.2 0.2 0.1], 'LineStyle', '-')
plot(x, [SingleRun1000_95pctup; SingleRun1000_95pctdown], '--r', 'LineWidth', 2)
plot(x, allCores_mixLog1000(:,2),  '-k', 'LineWidth',2)
plot(x, lognorm_BIGMACS.Var2, '--k', 'LineWidth',2)
title("1000 yr restriction")
xlim([0 6])
ylim(ylimits)
xlabel("nSR")
subplot(3,3,4)
hold on
plot(x, SR_MixLogNorm1Run0_high, 'Color', [1 0 0 0.1], 'LineStyle', '-')
plot(x, [SingleRun0_high_95pctup; SingleRun0_high_95pctdown], '--r', 'LineWidth', 2)
plot(x, highSRCores_mixLog0(:,2),  '-k', 'LineWidth',2)
plot(x, lognorm_BIGMACS.Var2, '--k', 'LineWidth',2)
title("0 yr restriction")
xlim([0 6])
ylim(ylimits)
xlabel("nSR")
subplot(3,3,5)
hold on
plot(x, SR_MixLogNorm1Run500_high, 'Color', [1 0 0 0.1], 'LineStyle', '-')
plot(x, [SingleRun500_high_95pctup; SingleRun500_high_95pctdown], '--r', 'LineWidth', 2)
plot(x, highSRCores_mixLog500(:,2),  '-k', 'LineWidth',2)
plot(x, lognorm_BIGMACS.Var2, '--k', 'LineWidth',2)
title("500 yr restriction")
xlim([0 6])
ylim(ylimits)
xlabel("nSR")
subplot(3,3,6)
hold on
plot(x, SR_MixLogNorm1Run1000_high, 'Color', [1 0 0 0.1], 'LineStyle', '-')
plot(x, [SingleRun1000_high_95pctup; SingleRun1000_high_95pctdown], '--r', 'LineWidth', 2)
plot(x, highSRCores_mixLog1000(:,2),  '-k', 'LineWidth',2)
plot(x, lognorm_BIGMACS.Var2, '--k', 'LineWidth',2)
title("1000 yr restriction")
xlim([0 6])
ylim(ylimits)
xlabel("nSR")
subplot(3,3,7)
hold on
plot(x, SR_MixLogNorm1Run0_low, 'Color', [0 0 1 0.1], 'LineStyle', '-')
plot(x, [SingleRun0_low_95pctup; SingleRun0_low_95pctdown], '--r', 'LineWidth', 2)
plot(x, lowSRCores_mixLog0(:,2),  '-k', 'LineWidth',2)
plot(x, lognorm_BIGMACS.Var2, '--k', 'LineWidth',2)
title("0 yr restriction")
xlim([0 6])
ylim(ylimits)
xlabel("nSR")
subplot(3,3,8)
hold on
plot(x, SR_MixLogNorm1Run500_low, 'Color', [0 0 1 0.1], 'LineStyle', '-')
plot(x, [SingleRun500_low_95pctup; SingleRun500_low_95pctdown], '--r', 'LineWidth', 2)
plot(x, lowSRCores_mixLog500(:,2),  '-k', 'LineWidth',2)
plot(x, lognorm_BIGMACS.Var2, '--k', 'LineWidth',2)
title("500 yr restriction")
xlim([0 6])
ylim(ylimits)
xlabel("nSR")
subplot(3,3,9)
hold on
plot(x, SR_MixLogNorm1Run1000_low, 'Color', [0 0 1 0.1], 'LineStyle', '-')
plot(x, [SingleRun1000_low_95pctup; SingleRun1000_low_95pctdown], '--r', 'LineWidth', 2)
plot(x, lowSRCores_mixLog1000(:,2),  '-k', 'LineWidth',2)
plot(x, lognorm_BIGMACS.Var2, '--k', 'LineWidth',2)
title("1000 yr restriction")
xlim([0 6])
ylim(ylimits)
xlabel("nSR")

figure;
subplot(3,1,1)
plot(x, lognorm_BIGMACS.Var2, '--k', 'LineWidth',2)
ylim([0 1])
xlim([0 6])
subplot(3,1,3)
plot(x, lowSRCores_mixLog500(:,2),  '-b', 'LineWidth',2)
hold on
plot(x, highSRCores_mixLog500(:,2),  '-r', 'LineWidth',2)
ylim([0 1])
xlim([0 6])
subplot(3,1,2)
plot(x, highSRCores_mixLog500(:,2),  '-r', 'LineWidth',2)
ylim([0 1])
xlim([0 6])

figure;
hold on
plot(x, highSRCores_mixLog500(:,2), 'r-', 'LineWidth', 2, "DisplayName", "SR > 8cm/kyr")
plot(x, lowSRCores_mixLog500(:,2),  '-b', 'LineWidth',2, "DisplayName", "SR < 8cm/kyr")
plot(x, lognorm_BIGMACS.Var2, '--k', 'LineWidth',2)
ylim([0 1])
xlim([0 6])
xlabel("nSR")
title("Mix Log Norm Comparison, 500yr Restriction")

%% GAMMA DISTRIBUTION FIGURES (INVERSE SR)
%Set up inverse nSR x-values
invx = sort(1./x);
invXbinEdges = [0:0.1:15];

%Fit gamma distributions to inverse nSR values with all cores
% [xgamma, gammaAllRunsSR, gammaAllRunsinvSR, alpha1, invSRhist1] = fitGamma2invSR(dataT.nSRcounts, allCoresLog, invXbinEdges);
% [gamma1RuninvSR,phat95up1,phat95low1,phats] = SingleRunGammas(dataT.nSRcounts, allCoresLog, 400, x, 1, 3, 1, 0);
% 
% [xgamma2, gammaAllRunsSR2, gammaAllRunsinvSR2, alpha2,invSRhist2] = fitGamma2invSR(dataT.nSRcounts500, allCoresLog, invXbinEdges);
% [gamma1RuninvSR2,phat95up2,phat95low2,phats2] = SingleRunGammas(dataT.nSRcounts500, allCoresLog, 400, x, 1, 3, 1, 0);

[xgamma2, gammaAllRunsSR3, gammaAllRunsinvSR3, alpha3,invSRhist3] = fitGamma2invSR(dataT.nSRcounts1000, allCoresLog, invXbinEdges);
[gamma1RuninvSR3,phat95up3,phat95low3,phats3] = SingleRunGammas(dataT.nSRcounts1000, allCoresLog, 400, x, 1, 3, 1, 0);

% [xgamma3, gammaAllRunsSR4, gammaAllRunsinvSR4, alpha4,invSRhist4] = fitGamma2invSR(dataT.nSRcounts, lowSRCoresLog, invXbinEdges);
% [gamma1RuninvSR4,phat95up4,phat95low4,phats4] = SingleRunGammas(dataT.nSRcounts, lowSRCoresLog, 400, x, 1, 3, 1, 0);
% 
% [xgamma3, gammaAllRunsSR5, gammaAllRunsinvSR5, alpha5,invSRhist5] = fitGamma2invSR(dataT.nSRcounts500, lowSRCoresLog, invXbinEdges);
% [gamma1RuninvSR5,phat95up5,phat95low5,phats5] = SingleRunGammas(dataT.nSRcounts500, lowSRCoresLog, 400, x, 1, 3, 1, 0);

[xgamma3, gammaAllRunsSR6, gammaAllRunsinvSR6, alpha6,invSRhist6] = fitGamma2invSR(dataT.nSRcounts1000, lowSRCoresLog, invXbinEdges);
[gamma1RuninvSR6,phat95up6,phat95low6,phats6] = SingleRunGammas(dataT.nSRcounts1000, lowSRCoresLog, 400, x, 1, 3, 1, 0);

% [xgamma3, gammaAllRunsSR7, gammaAllRunsinvSR7, alpha7,invSRhist7] = fitGamma2invSR(dataT.nSRcounts, highSRCoresLog, invXbinEdges);
% [gamma1RuninvSR7,phat95up7,phat95low7,phats7] = SingleRunGammas(dataT.nSRcounts, highSRCoresLog, 400, x, 1, 3, 1, 0);
% 
% [xgamma3, gammaAllRunsSR8, gammaAllRunsinvSR8, alpha8,invSRhist8] = fitGamma2invSR(dataT.nSRcounts500, highSRCoresLog, invXbinEdges);
% [gamma1RuninvSR8,phat95up8,phat95low8,phats8] = SingleRunGammas(dataT.nSRcounts500, highSRCoresLog, 400, x, 1, 3, 1, 0);

[xgamma3, gammaAllRunsSR9, gammaAllRunsinvSR9, alpha9,invSRhist9] = fitGamma2invSR(dataT.nSRcounts1000, highSRCoresLog, invXbinEdges);
[gamma1RuninvSR9,phat95up9,phat95low9,phats9] = SingleRunGammas(dataT.nSRcounts1000, highSRCoresLog, 400, x, 1, 3, 1, 0);

%% Plot impact of sampling resolution on pdf

figure;
subplot(3,2,1)
histogram("BinCounts", invSRhist1, "BinEdges", invXbinEdges, "Normalization", "pdf", "FaceColor",'k', "FaceAlpha",0.2)
xlabel("inverse nSR")
xlim([0 6])
subplot(3,2,2)
histogram("BinCounts", allCores0_diffData, "BinEdges", agediffsBinEdges, "FaceColor",'k', "FaceAlpha",0.2)
xlabel("Age Diff (yr)")
ylim([0 4000])

subplot(3,2,3)
histogram("BinCounts", invSRhist2, "BinEdges", invXbinEdges, "Normalization", "pdf", "FaceColor",'k', "FaceAlpha",0.2)
xlabel("inverse nSR")
xlim([0 6])
subplot(3,2,4)
histogram("BinCounts", allCores500_diffData, "BinEdges", agediffsBinEdges, "FaceColor",'k', "FaceAlpha",0.2)
xlabel("Age Diff (yr)")
ylim([0 4000])

subplot(3,2,5)
histogram("BinCounts", invSRhist2, "BinEdges", invXbinEdges, "Normalization", "pdf", "FaceColor",'k', "FaceAlpha",0.2)
xlabel("inverse nSR")
xlim([0 6])
subplot(3,2,6)
histogram("BinCounts", allCores1000_diffData, "BinEdges", agediffsBinEdges, "FaceColor",'k', "FaceAlpha",0.2)
xlabel("Age Diff (yr)")
ylim([0 4000])

%% Plot histogram of allcores 500res data in invSR with inset of resolution histogram

hist_with_inset0= figure;
histogram("BinCounts", invSRhist1, "BinEdges", invXbinEdges, "FaceColor",'k', "FaceAlpha",0.2)
xlim([0 5])
xlabel("Inverse nSR")
ylabel("Sediment Length (cm)")
ylim([0 1800])
handaxes2 = axes("Position", [0.52 0.42 0.3 0.3]);
histogram("BinCounts", allCores0_diffData, "BinEdges", agediffsBinEdges, "FaceColor",'k', "FaceAlpha",0.2)
xlabel("Age Diff (yr)")
ylim([0 3100])
fontsize(hist_with_inset0, 20, 'points')

hist_with_inset= figure;
histogram("BinCounts", invSRhist2, "BinEdges", invXbinEdges, "FaceColor",'k', "FaceAlpha",0.2)
xlim([0 5])
xlabel("Inverse nSR")
ylabel("Sediment Length (cm)")
handaxes3 = axes("Position", [0.52 0.42 0.3 0.3]);
histogram("BinCounts", allCores500_diffData, "BinEdges", agediffsBinEdges, "FaceColor",'k', "FaceAlpha",0.2)
xlabel("Age Diff (yr)")
ylim([0 3100])
fontsize(hist_with_inset, 20, 'points')

%% Plot results
figure;
subplot(3,3,1)
hold on
plot(invx, gamma1RuninvSR, 'Color', [0.2 0.2 0.2 0.1], 'LineStyle', '-')
plot(invx, gammaAllRunsinvSR, '-k', 'LineWidth', 2)
title("0 year restriction")
 text(3, 0.8, {"allRuns = "+ num2str(alpha1, '%.2f'),"95% low = " + num2str(phat95low1,'%.2f'),"95% up = " + num2str(phat95up1,'%.2f')})
ylim([0 1])
xlim([0 5])
xlabel("Inverse nSR")
subplot(3,3,2)
hold on
plot(invx, gamma1RuninvSR2, 'Color', [0.2 0.2 0.2 0.1], 'LineStyle', '-')
plot(invx, gammaAllRunsinvSR2, '-k', 'LineWidth', 2)
title("500 year restriction")
 text(3, 0.8, {"allRuns = "+ num2str(alpha2, '%.2f'),"95% low = " + num2str(phat95low2,'%.2f'),"95% up = " + num2str(phat95up2,'%.2f')})
ylim([0 1])
xlim([0 5])
xlabel("Inverse nSR")
subplot(3,3,3)
hold on
plot(invx, gamma1RuninvSR3, 'Color', [0.2 0.2 0.2 0.1], 'LineStyle', '-')
plot(invx, gammaAllRunsinvSR3, '-k', 'LineWidth', 2)
title("1000 year restriction")
 text(3, 0.8, {"allRuns = "+ num2str(alpha3, '%.2f'),"95% low = " + num2str(phat95low3,'%.2f'),"95% up = " + num2str(phat95up3,'%.2f')})
ylim([0 1])
xlim([0 5])
xlabel("Inverse nSR")

subplot(3,3,4)
hold on
plot(invx, gamma1RuninvSR7, 'Color', [0.2 0.2 0.2 0.1], 'LineStyle', '-')
plot(invx, gammaAllRunsinvSR7, '-r', 'LineWidth', 2)
 text(3, 0.8, {"allRuns = "+ num2str(alpha7, '%.2f'),"95% low = " + num2str(phat95low7,'%.2f'),"95% up = " + num2str(phat95up7,'%.2f')})
ylim([0 1])
xlim([0 5])
xlabel("Inverse nSR")
subplot(3,3,5)
hold on
plot(invx, gamma1RuninvSR8, 'Color', [0.2 0.2 0.2 0.1], 'LineStyle', '-')
plot(invx, gammaAllRunsinvSR8, '-r', 'LineWidth', 2)
 text(3, 0.8, {"allRuns = "+ num2str(alpha8, '%.2f'),"95% low = " + num2str(phat95low8,'%.2f'),"95% up = " + num2str(phat95up8,'%.2f')})
ylim([0 1])
xlim([0 5])
xlabel("Inverse nSR")
subplot(3,3,6)
hold on
plot(invx, gamma1RuninvSR9, 'Color', [0.2 0.2 0.2 0.1], 'LineStyle', '-')
plot(invx, gammaAllRunsinvSR9, '-r', 'LineWidth', 2)
 text(3, 0.8, {"allRuns = "+ num2str(alpha9, '%.2f'),"95% low = " + num2str(phat95low9,'%.2f'),"95% up = " + num2str(phat95up9,'%.2f')})
ylim([0 1])
xlim([0 5])
xlabel("Inverse nSR")

subplot(3,3,7)
hold on
plot(invx, gamma1RuninvSR4, 'Color', [0.2 0.2 0.2 0.1], 'LineStyle', '-')
plot(invx, gammaAllRunsinvSR4, '-b', 'LineWidth', 2)
 text(3, 0.8, {"allRuns = "+ num2str(alpha4, '%.2f'),"95% low = " + num2str(phat95low4,'%.2f'),"95% up = " + num2str(phat95up4,'%.2f')})
ylim([0 1])
xlim([0 5])
xlabel("Inverse nSR")
subplot(3,3,8)
hold on
plot(invx, gamma1RuninvSR5, 'Color', [0.2 0.2 0.2 0.1], 'LineStyle', '-')
plot(invx, gammaAllRunsinvSR5, '-b', 'LineWidth', 2)
 text(3, 0.8, {"allRuns = "+ num2str(alpha5, '%.2f'),"95% low = " + num2str(phat95low5,'%.2f'),"95% up = " + num2str(phat95up5,'%.2f')})
ylim([0 1])
xlim([0 5])
xlabel("Inverse nSR")
subplot(3,3,9)
hold on
plot(invx, gamma1RuninvSR6, 'Color', [0.2 0.2 0.2 0.1], 'LineStyle', '-')
plot(invx, gammaAllRunsinvSR6, '-b', 'LineWidth', 2)
 text(3, 0.8, {"allRuns = "+ num2str(alpha6, '%.2f'),"95% low = " + num2str(phat95low6,'%.2f'),"95% up = " + num2str(phat95up6,'%.2f')})
ylim([0 1])
xlim([0 5])
xlabel("Inverse nSR")

%%  Compare Gamma distributions from high and low SR

%Plot results of full sampling fitting
figure;
hold on
plot(invx, gammaAllRunsinvSR6, '-b', 'LineWidth', 2, "DisplayName", "SR < 8cm/kyr")
plot(invx, gammaAllRunsinvSR9, '-r', 'LineWidth', 2, "DisplayName", "SR > 8cm/kyr")
xlim([0 5])
xlabel("Inverse nSR")
legend()

%Plot results of single run fitting 
figure;
hold on
plot(invx, gamma1RuninvSR6, 'Color', [0 0 1 0.1], 'LineStyle', '-')
plot(invx, gamma1RuninvSR9, 'Color', [1 0 0 0.1], 'LineStyle', '-')
 text(3, 0.8, {"SR > 8cm/kyr","95% CI = [" + num2str(phat95low9,'%.2f') + ", " + num2str(phat95up9,'%.2f') + "]","", "SR < 8cm/kyr", " 95% CI = [" + num2str(phat95low6,'%.2f') + ", " + num2str(phat95up6,'%.2f') + "]"})
xlim([0 5])
xlabel("Inverse nSR")

%% Plot effects of data resolution on highSR dataset
figure;
hold on
plot(invx, gammaAllRunsinvSR9, 'g', 'LineWidth', 2, "DisplayName", "1000yr; alpha = " + num2str(alpha9, '%.2f'))
plot(invx, gammaAllRunsinvSR8, 'r', 'LineWidth', 2, "DisplayName", "500yr; alpha = " + num2str(alpha8, '%.2f'))
plot(invx, gammaAllRunsinvSR7, 'k', 'LineWidth', 2, "DisplayName", "0yr; alpha = " + num2str(alpha7, '%.2f'))
xlim([0 5])
xlabel("Inverse nSR")
title("HighSR")
legend()

figure;
hold on
plot(invx, gamma1RuninvSR9, 'Color', [0 1 0 0.1], 'LineStyle', '-')
plot(invx, gamma1RuninvSR8, 'Color', [1 0 0 0.1], 'LineStyle', '-')
plot(invx, gamma1RuninvSR7, 'Color', [0 0 0 0.1], 'LineStyle', '-')
xlim([0 5])
%text(3, 0.8, {"95% low = " + num2str(phat95low9,'%.2f'),"95% up = " + num2str(phat95up9,'%.2f'), "", "95% low = " + num2str(phat95low8,'%.2f'),"95% up = " + num2str(phat95up8,'%.2f'), "", "95% low = " + num2str(phat95low7,'%.2f'),"95% up = " + num2str(phat95up7,'%.2f')})
text(3, 0.8, {"1000 year; ", "95% CI = [" + num2str(phat95low9,'%.2f')+ ", " + num2str(phat95up9,'%.2f')+ "]", "", "500 year;", "95% CI = [" + num2str(phat95low8,'%.2f') + ", " + num2str(phat95up8,'%.2f') + "]", "", "0 year",  "95% CI = [" + num2str(phat95low7,'%.2f') + ", " + num2str(phat95up7,'%.2f') + "]"})
xlabel("Inverse nSR")
title("HighSR")

figure;
hold on
plot(invx, gammaAllRunsinvSR6, 'g', 'LineWidth', 2, "DisplayName", "1000yr; alpha = " + num2str(alpha6, '%.2f'))
plot(invx, gammaAllRunsinvSR5, 'b', 'LineWidth', 2, "DisplayName", "500yr; alpha = " + num2str(alpha5, '%.2f'))
plot(invx, gammaAllRunsinvSR4, 'k', 'LineWidth', 2, "DisplayName", "0yr; alpha = " + num2str(alpha4, '%.2f'))
xlim([0 5])
xlabel("Inverse nSR")
title("LowSR")
legend()

figure;
hold on
plot(invx, gamma1RuninvSR6, 'Color', [0 1 0 0.1], 'LineStyle', '-')
plot(invx, gamma1RuninvSR5, 'Color', [0 0 1 0.1], 'LineStyle', '-')
plot(invx, gamma1RuninvSR4, 'Color', [0 0 0 0.1], 'LineStyle', '-')
xlim([0 5])
text(3, 0.8, {"95% low = " + num2str(phat95low6,'%.2f'),"95% up = " + num2str(phat95up6,'%.2f'), "", "95% low = " + num2str(phat95low5,'%.2f'),"95% up = " + num2str(phat95up5,'%.2f'), "", "95% low = " + num2str(phat95low4,'%.2f'),"95% up = " + num2str(phat95up4,'%.2f')})
xlabel("Inverse nSR")
title("LowSR")

%% Plot gamma vs mix log norm in nSR space

figure;
hold on
plot(xgamma3, gammaAllRunsSR7, 'k--', 'LineWidth', 2, "DisplayName", "gamma0")
plot(xgamma3, gammaAllRunsSR8, 'r--', 'LineWidth',2, "DisplayName", "gamma500")
plot(xgamma3, gammaAllRunsSR9, 'g--', 'LineWidth',2, "DisplayName", "gamma1000")
plot(x, highSRCores_mixLog500(:,2),  '-r', 'LineWidth',2, "DisplayName", "MixLogNorm500")
xlim([0 5])             
xlabel("nSR")
title("HighSR 500")
legend()

%% compare gamma and mix log to histogram
figure;
hold on
plot(xgamma3, gammaAllRunsSR9, 'r--', 'LineWidth',2, "DisplayName", "Gamma - This Study")
plot(x, highSRCores_mixLog1000(:,2),  '-r', 'LineWidth',2, "DisplayName", "MixLogNorm1000")
plot(lognorm_BIGMACS.Var1, lognorm_BIGMACS.Var2,  'k-', 'LineWidth',2, "DisplayName", "BIGMACS")
histogram(highSRCores1000_histData, "Normalization", "pdf", "FaceColor",'k', "FaceAlpha",0.2, "DisplayName", "High SR Cores Data")
xlim([0 5])             
xlabel("nSR")
title("HighSR 1000")
legend()

%% 

%% compare gamma and mix log norm to histogram in invNSR space

[invx2, MLNinverse_high1000] = invSRtoSR( x, highSRCores_mixLog1000(:,2));
[invxBM, BMinverse] = invSRtoSR(lognorm_BIGMACS.Var1, lognorm_BIGMACS.Var2);

figure; 
hold on
plot(invx, gammaAllRunsinvSR9, 'r--', 'LineWidth', 2, "DisplayName", "500yr; alpha = " + num2str(alpha9, '%.2f'))
%plot(invx, gamma1RuninvSR8, 'Color', [1 0 0 0.1], 'LineStyle', '-')
plot(invx2, MLNinverse_high1000, 'r-', 'LineWidth', 2, "DisplayName", "MixLogNorm")
plot(invxBM, BMinverse, 'k-','LineWidth', 2, "DisplayName", "BIGMACS")
histogram("BinCounts", invSRhist9, "BinEdges", invXbinEdges, "Normalization","pdf", "FaceColor",'k', "FaceAlpha",0.2, "DisplayName", "High SR Cores Data")
xlim([0 5])
xlabel("inverse nSR")
title("HighSR 1000")
legend()