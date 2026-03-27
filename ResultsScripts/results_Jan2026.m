%% Add necessary paths
addpath("../Functions/")

%Set text size
set(groot, 'DefaultAxesFontSize', 12)
set(groot, 'DefaultTextFontSize', 12)

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
[invGam_BIGMACS, ~, BM.invgamfit] = fitInvGamma(BM.hist, BM.nSR.x, 343);
[BM.lnSR.x, BM.lnSR.px] = px_to_pfx(BM.nSR.x, BM.nSR.px, @log);
BMpdf.numParams = 6;
BMpdf.x = BM.lnSR.x;
BMpdf.px = BM.lnSR.px;

%% Load data and fits
%load("../Results/dataT_All1_RLGtrue_DS0p05_Dec9_fit28Feb26_depthweight_400R");
load("../Results/dataT_All1_RLGtrue_DS0p05_Dec9_fit26Mar26_noweight_SR_400R.mat")

%% Take data of interest
%----- Find corenames, lats, longs, depths of cores included
dataTbl = d.dataT(d.S1.chooseLog, :);
%----- Make map with locations denoted as red squares

% Define tick locations
lonTicks = -180:20:180;
latTicks = -40:20:40;
tickLen = 2; % in degrees, adjust to taste


figure;

% Map panel
subplot(2,2,[1 2])
ax = axesm('mercator', 'MapLatLimit', [-45 45], 'MapLonLimit', [-180 180], ...
    'MeridianLabel', 'on', 'ParallelLabel', 'on', 'Frame', 'off'...
    ,'LabelFormat', 'none', 'MLabelLocation', 20, 'PLabelLocation', 20,...
    'MLabelParallel', -45);
hold on
geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', [0.8 0.8 0.8])
load coastlines coastlat coastlon
plotm(coastlat, coastlon, 'Color', 'k')

tightmap

plotm(dataTbl.lats, dataTbl.longs, 'ks', 'LineWidth', 1)
set(gca, 'OuterPosition', [0 0.4 1 0.6])  % [left, bottom, width, height]
text(0.02, 0.9, 'A', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')


% Draw on ticks
xl = xlim;
yl = ylim;
xTickLen = 0.05;
yTickLen = 0.05;

proj = gcm;

[lonTick_x, lonTick_y] = projfwd(proj, repmat(mapLatLim(1), size(lonTicks)), lonTicks);
[latTick_x, latTick_y] = projfwd(proj, latTicks, repmat(mapLonLim(1), size(latTicks)));

% Draw longitude ticks at bottom and top
for i = 1:length(lonTicks)
    line([lonTick_x(i) lonTick_x(i)], [yl(1) yl(1)+yTickLen], 'Color', 'k', 'LineWidth', 0.5)
    line([lonTick_x(i) lonTick_x(i)], [yl(2)-yTickLen yl(2)], 'Color', 'k', 'LineWidth', 0.5)
end

% Draw latitude ticks at left and right
for i = 1:length(latTicks)
    line([xl(1) xl(1)+xTickLen], [latTick_y(i) latTick_y(i)], 'Color', 'k', 'LineWidth', 0.5)
    line([xl(2)-xTickLen xl(2)], [latTick_y(i) latTick_y(i)], 'Color', 'k', 'LineWidth', 0.5)
end

%Rotate longitude labels
xticklabelHandles = findall(gcf,'Type','text','Tag','MLabel');
nLabeles = numel(xticklabelHandles);
xtickPosition = cell2mat(get(xticklabelHandles,'Position'));
hold(ax,'on')
xTickLabels = get(xticklabelHandles, 'String'); 
xTickLabels = cellfun(@(c)c(2),xTickLabels);
set(xticklabelHandles, {'String'}, xTickLabels)
% Rotate text
set(xticklabelHandles,'Rotation', 90)
set(xticklabelHandles, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'Middle')
% Shift labels downward
xtickPosition(:,2) = xtickPosition(:,2) - range(ax.YLim)*.03;
set(xticklabelHandles, {'Position'}, mat2cell(xtickPosition, ones(nLabeles,1),3))
delete(anchors)
tightmap

%Manipulate y tick labels
yticklabelHandles = findall(gcf,'Type','text','Tag','PLabel');
nyLabeles = numel(yticklabelHandles);
ytickPosition = cell2mat(get(yticklabelHandles,'Position'));
ytickPosition(:,1) = ytickPosition(:,1) + range(ax.XLim)*.02;
set(yticklabelHandles, {'Position'}, mat2cell(ytickPosition, ones(nyLabeles,1),3))
tightmap

%%

subplot(2,2,3)
histogram(dataTbl.depths./1000, 0:0.5:6, 'FaceColor', [0.8 0.8 0.8])
xlim([0 6])
xlabel("Depth (km)")
ylabel("Counts")
xticks(0:1:6)
ylims = ylim; ylim(ylims.*1.1)
set(gca, 'Position', [0.13 0.25 0.33 0.25])
text(0.02, 0.9, 'B', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

subplot(2,2,4)
histogram(dataTbl.meanSR, 0:5:90, 'FaceColor', [0.8 0.8 0.8])
xlabel('Average SR (cm/kyr)')
ylabel('Counts')
xlim([0 90])
ylims = ylim; ylim(ylims.*1.1)
set(gca, 'Position', [0.57 0.25 0.33 0.25])
text(0.02, 0.9, 'C', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')



%% Plot age modes 
%Find out which cores passed filtering
usedLog = ~isnan(d.dataT.meanSR);

plotAgeModes(d.S1.chooseLog,d.S1.chooseLog, d.dataT.ageModes, d.dataT.cores)

%% %Create table of useful information
dStrus = {d.S1.BMedian, d.S1.BChIR, d.S1.RSR0IR, d.S1.RSR500IR, d.S1.RSR1000IR};
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

%% Plot of BIGMACS and BMedian approaches with normalized BIGMACS histogram
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
BMedian_logHC = histcounts(log(d.S1.BMedian.weightedC), logBinEdges)./d.S1.fitS.DeterministicRun.weightInflator;
normF = 1./sum((BMedian_logHC.*uniquetol(diff(logBinEdges), 1e-5)));
histogram('BinCounts', BMedian_logHC,'BinEdges', logBinEdges, 'FaceColor', 'r', 'FaceAlpha', 0.1);
hold on
plot(BM.lnSR.x, BM.lnSR.px*(1/normF), 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--')
box on
%xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BMedian")
xlabel("log(NSR)")

%saveas(fA, "../Figures/BMvsBMedian_normalizedPDF.jpg")

%% Plot of BIGMACS and BMedian with normalized pdfs and histogram of dt

load("../Results/Lin2014_dts.mat")

%Get weighted histogram of dts for BIGMACS
BM_dt_weightedBC = makeWeightedBinCounts(asr(:,2), asr(:,3), 0:0.1:10);

nsubs = 4;
fA = figure;
fA.Position = [100, 100, 1000, 400]; % Wide and short: [left, bottom, width, height]

subplot(2,2,1)
BM.logHC = histcounts(log(BM.hist), logBinEdges);
normF = 1./sum((BM.logHC.*uniquetol(diff(logBinEdges), 1e-5)));
hold on
histogram('BinCounts', BM.logHC, 'BinEdges', logBinEdges, 'FaceColor', [0.8 0.8 0.8])
plot(BM.lnSR.x, BM.lnSR.px.*(1/normF), 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--')
text(0.02, 0.9, 'A', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')
box on
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
xlabel("log(NSR)")
subplot(2,2,2)
hold on
histogram('BinCounts', BM_dt_weightedBC, 'BinEdges', 0:0.1:10,  'FaceColor', [0.8 0.8 0.8])
text(0.02, 0.9, 'B', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')
box on
ylabel("cm")
xlabel("\Deltat (kyr)")
xlim([0 5])
ylim([0 2000])
subplot(2,2,3)
BMedian_logHC = histcounts(log(d.S1.BMedian.weightedC), logBinEdges);
normF = 1./sum((BMedian_logHC.*uniquetol(diff(logBinEdges), 1e-5)));
histogram('BinCounts', BMedian_logHC,'BinEdges', logBinEdges,  'FaceColor', [0.8 0.8 0.8]);
hold on
plot(BM.lnSR.x, BM.lnSR.px*(1/normF), 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--')
text(0.02, 0.9, 'C', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')
box on
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
xlabel("log(NSR)")
subplot(2, 2, 4)
histogram('BinCounts', d.S1.BMedian.agediffsWC, 'BinEdges', 0:0.1:10,  'FaceColor', [0.8 0.8 0.8])
text(0.02, 0.9, 'D', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')
box on
ylabel("cm")
xlabel("\Deltat (kyr)")
xlim([0 5])
ylim([0 2000])

%% Plot of pooled samples + fits

AR = d.S1.RSR500AR;
IR = d.S1.RSR500IR;

RGB = orderedcolors("gem");
H = rgb2hex(RGB);

figure;
subplot(2,1,1)
hold on
AR_logHC = histcounts(log(AR.weightedC), logBinEdges);
normF = 1./sum((AR_logHC.*uniquetol(diff(logBinEdges), 1e-5)));
histogram('BinCounts', AR_logHC, 'BinEdges', logBinEdges, 'HandleVisibility', 'off', 'FaceColor', '[0.8 0.8 0.8]')
plot(AR.LN.lnSR.x, AR.LN.lnSR.px*(1/normF), 'LineWidth', 1, 'DisplayName','LN', 'Color', 'r')
plot(AR.MLN.lnSR.x, AR.MLN.lnSR.px*(1/normF), 'LineWidth', 1, 'DisplayName','MLN', 'Color', H(5))
plot(AR.Gam.lnSR.x, AR.Gam.lnSR.px*(1/normF),  'LineWidth', 1, 'DisplayName','G', 'Color', H(4) )
plot(AR.invGam.lnSR.x, AR.invGam.lnSR.px*(1/normF),  'LineWidth', 1, 'DisplayName','IG', 'Color', 'b')
plot(BM.lnSR.x, BM.lnSR.px.*(1/normF), 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'DisplayName', 'BIGMACS')
xlim([-2.5 2.5])
xlabel("log(NSR)")
ylabel("Pooled Counts")
legend()
ylim([0 1.2e5])
ylims = ylim;
text(0.02, 0.9, 'A', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')


%And the individual runs as well
subplot(2,1,2)
hold on;
AR_logHC = histcounts(log(AR.weightedC), logBinEdges);
normF = 1./sum((AR_logHC.*uniquetol(diff(logBinEdges), 1e-5)));
histogram('BinCounts', AR_logHC, 'BinEdges', logBinEdges, 'HandleVisibility', 'off', 'FaceColor', '[0.8 0.8 0.8]')
plot(IR.LN.lnSR.x, IR.LN.lnSR.px(:,1).*(1/normF), 'Color', [0 0 0 0.1], 'DisplayName', 'Individual Run LN')
plot(IR.LN.lnSR.x, IR.LN.lnSR.px(:,2:end).*(1/normF), 'Color', [0 0 0 0.1], 'HandleVisibility', 'off')
plot(AR.LN.lnSR.x, AR.LN.lnSR.px.*(1/normF), 'Color', "r", 'LineWidth', 1, 'DisplayName', 'LN')
%plot(BM.lnSR.x, BM.lnSR.px.*(1/normF), 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'DisplayName', 'BIGMACS')
xlim([-2.5 2.5])
ylim(ylims)
xlabel("log(NSR)")
ylabel("Pooled Counts")
legend()
text(0.02, 0.9, 'B', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

%% Plot of pooled samples + fits on inverse NSR axis

AR = d.S1.RSR500AR;
IR = d.S1.RSR500IR;
invBinEdges = 0:0.1:6;

f_inv = @(x)(1./x);
[AR.LN.invSR.x, AR.LN.invSR.px] = px_to_pfx(AR.LN.nSR.x, AR.LN.nSR.px, f_inv);
[AR.invGam.invSR.x, AR.invGam.invSR.px] = px_to_pfx(AR.invGam.nSR.x, AR.invGam.nSR.px, f_inv);


figure;
hold on
AR_invHC = histcounts(1./AR.weightedC, invBinEdges);
normF = 1./sum((AR_invHC.*uniquetol(diff(invBinEdges), 1e-5)));
histogram('BinCounts', AR_invHC, 'BinEdges', invBinEdges, 'HandleVisibility', 'off', 'FaceColor', '[0.8 0.8 0.8]')
plot(AR.LN.invSR.x, AR.LN.invSR.px*(1/normF), 'r', 'LineWidth', 1, 'DisplayName','LN')
plot(AR.invGam.invSR.x, AR.invGam.invSR.px*(1/normF),  'LineWidth', 1, 'DisplayName','IG', 'Color', 'b')
xlim([0 6])
xlabel("inverse NSR")
ylabel("Pooled Counts")
legend()
%ylim([0 1.2e5])
ylims = ylim;

%% Double check the pooled weighted histograms match the pooled weighted histograms of individual runs

% IR_WCpool = [];
% for i = 1:length(IR.weightedC)
%     IR_WCpool = [IR_WCpool, IR.weightedC{i}];
% end
% 
% figure;
% subplot(2,1,1)
% histogram('BinCounts', AR_logHC, 'BinEdges', logBinEdges, 'HandleVisibility', 'off', 'FaceColor', '[0.8 0.8 0.8]')
% subplot(2,1,2)
% IR_logHC = histcounts(log(IR_WCpool), logBinEdges);
% histogram('BinCounts', IR_logHC, 'BinEdges', logBinEdges, 'HandleVisibility', 'off', 'FaceColor', '[0.8 0.8 0.8]')
% 


%% All sampling approaches histogram of median with 68th percentile bars + dt histogram
nsubs = 4;
numruns = length(d.S1.BChIR.OneRunDatas);
fB = figure;
fB.Position = [100, 100, 1000, 800]; % Wide and short: [left, bottom, width, height]
subplot(4,2,1)
A = sort(d.S1.BChIR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceColor', '[0.8 0.8 0.8]')
hold on
%plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
xlabel("log(NSR)")
ylabel("cm")
text(0.02, 0.9, 'A', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

subplot(4,2,2)
B = sort(d.S1.BChIR.agediffsWbC', 1);
hold on
box on
histogram('BinCounts', B(numruns*0.5, :), 'BinEdges', 0:0.1:10, 'FaceColor', '[0.8 0.8 0.8]')
errorbar(0.05:0.1:9.95, B(numruns*0.5, :),(B(numruns*0.5, :)-B(numruns*0.16, :)),(B(numruns*0.84, :)-B(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
xlim([0 6])
ylim([0 2000])
xlabel("\Deltat (kyr)")
ylabel("cm")
text(0.02, 0.9, 'B', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')


subplot(4,2,3)
A = sort(d.S1.RSR0IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceColor', '[0.8 0.8 0.8]')
hold on
%plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])

ylabel("cm")
xlabel("log(nSR)")
text(0.02, 0.9, 'C', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

subplot(4,2,4)
B = sort(d.S1.RSR0IR.agediffsWbC', 1);
hold on
box on
histogram('BinCounts', B(numruns*0.5, :), 'BinEdges', 0:0.1:10, 'FaceColor', '[0.8 0.8 0.8]')
errorbar(0.05:0.1:9.95, B(numruns*0.5, :),(B(numruns*0.5, :)-B(numruns*0.16, :)),(B(numruns*0.84, :)-B(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
xlim([0 6])
ylim([0 2000])

ylabel("cm")
xlabel("\Deltat (kyr)")
text(0.02, 0.9, 'D', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

subplot(4,2,5)
A = sort(d.S1.RSR500IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceColor', '[0.8 0.8 0.8]')
hold on
%plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])

ylabel("cm")
xlabel("log(NSR)")
text(0.02, 0.9, 'E', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

subplot(4,2,6)
B = sort(d.S1.RSR500IR.agediffsWbC', 1);
hold on
box on
histogram('BinCounts', B(numruns*0.5, :), 'BinEdges', 0:0.1:10, 'FaceColor', '[0.8 0.8 0.8]')
errorbar(0.05:0.1:9.95, B(numruns*0.5, :),(B(numruns*0.5, :)-B(numruns*0.16, :)),(B(numruns*0.84, :)-B(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
xlim([0 6])
ylim([0 2000])

ylabel("cm")
xlabel("\Deltat (kyr)")
text(0.02, 0.9, 'F', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

subplot(4,2,7)
A = sort(d.S1.RSR1000IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceColor', '[0.8 0.8 0.8]')
hold on
%plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
xlabel("log(NSR)")
ylabel("cm")
text(0.02, 0.9, 'G', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')


subplot(4,2,8)
B = sort(d.S1.RSR1000IR.agediffsWbC', 1);
hold on
box on
histogram('BinCounts', B(numruns*0.5, :), 'BinEdges', 0:0.1:10, 'FaceColor', '[0.8 0.8 0.8]')
errorbar(0.05:0.1:9.95, B(numruns*0.5, :),(B(numruns*0.5, :)-B(numruns*0.16, :)),(B(numruns*0.84, :)-B(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
xlim([0 6])
ylim([0 2000])
xlabel("\Deltat (kyr)")
ylabel("cm")
text(0.02, 0.9, 'H', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')


%% All sampling approaches histogram pooled
nsubs = 4;
numruns = length(d.S1.BChIR.OneRunDatas);
figure;
subplot(4,2,1)
A = sort(d.S1.BChIR.lnSRHistCounts, 1);
box on
histogram('BinCounts', sum(A,1), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
xlim([-2.5 2.5])
%ylim([0 3500])
xlabel("log(NSR)")
ylabel("cm")
title("BSamp [500-4000yr]")
subplot(4,2,2)
B = sort(d.S1.BChIR.agediffsWbC', 1);
hold on
box on
histogram('BinCounts', sum(B,1), 'BinEdges', 0:0.1:10, 'FaceAlpha', 0.1)
xlim([0 6])
%ylim([0 2000])
xlabel("\Deltat (kyr)")
ylabel("cm")
subplot(4,2,3)
A = sort(d.S1.RSR0IR.lnSRHistCounts, 1);
box on
histogram('BinCounts', sum(A,1), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
xlim([-2.5 2.5])
%ylim([0 3500])
xlabel("log(NSR)")
ylabel("cm")
%xlabel("log(nSR)")
title("RSR0")
subplot(4,2,4)
B = sort(d.S1.RSR0IR.agediffsWbC', 1);
hold on
box on
histogram('BinCounts', sum(B,1), 'BinEdges', 0:0.1:10, 'FaceAlpha', 0.1)
xlim([0 6])
%ylim([0 2000])
xlabel("\Deltat (kyr)")
ylabel("cm")


subplot(4,2,5)
A = sort(d.S1.RSR500IR.lnSRHistCounts, 1);
box on
histogram('BinCounts', sum(A,1), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
xlim([-2.5 2.5])
%ylim([0 3500])
xlabel("log(NSR)")
ylabel("cm")
title("RSR500")
subplot(4,2,6)
B = sort(d.S1.RSR500IR.agediffsWbC', 1);
hold on
box on
histogram('BinCounts', sum(B,1), 'BinEdges', 0:0.1:10, 'FaceAlpha', 0.1)
xlim([0 6])
%ylim([0 2000])
xlabel("\Deltat (kyr)")
ylabel("cm")

subplot(4,2,7)
A = sort(d.S1.RSR1000IR.lnSRHistCounts, 1);
box on
histogram('BinCounts', sum(A,1), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
xlim([-2.5 2.5])
%ylim([0 3500])
xlabel("log(NSR)")
ylabel("cm")
title("RSR1000")
subplot(4,2,8)
B = sort(d.S1.RSR1000IR.agediffsWbC', 1);
hold on
box on
histogram('BinCounts', sum(B,1), 'BinEdges', 0:0.1:10, 'FaceAlpha', 0.1)
xlim([0 6])
%ylim([0 2000])
xlabel("\Deltat (kyr)")
ylabel("cm")



%% BIC (taeheefix) of each method
%This is the BIC (with the correction to the likelihood) for each 

BIC_BMedian = round([d.S1.BMedian.LN.fitInfo.BICtaeheefix, d.S1.BMedian.MLN.fitInfo.BICtaeheefix, d.S1.BMedian.Gam.fitInfo.BICtaeheefix, d.S1.BMedian.invGam.fitInfo.BICtaeheefix])
BIC_BChAR = round([d.S1.BChAR.LN.fitInfo.BICtaeheefix, d.S1.BChAR.MLN.fitInfo.BICtaeheefix, d.S1.BChAR.Gam.fitInfo.BICtaeheefix, d.S1.BChAR.invGam.fitInfo.BICtaeheefix])
BIC_RSR0AR = round([d.S1.RSR0AR.LN.fitInfo.BICtaeheefix, d.S1.RSR0AR.MLN.fitInfo.BICtaeheefix, d.S1.RSR0AR.Gam.fitInfo.BICtaeheefix, d.S1.RSR0AR.invGam.fitInfo.BICtaeheefix])
BIC_RSR500AR = round([d.S1.RSR500AR.LN.fitInfo.BICtaeheefix, d.S1.RSR500AR.MLN.fitInfo.BICtaeheefix, d.S1.RSR500AR.Gam.fitInfo.BICtaeheefix, d.S1.RSR500AR.invGam.fitInfo.BICtaeheefix])
BIC_RSR1000AR = round([d.S1.RSR1000AR.LN.fitInfo.BICtaeheefix, d.S1.RSR1000AR.MLN.fitInfo.BICtaeheefix, d.S1.RSR1000AR.Gam.fitInfo.BICtaeheefix, d.S1.RSR1000AR.invGam.fitInfo.BICtaeheefix])

%% Standard Deviation of NSR counts for each approach

std(d.S1.BMedian.weightedC)
std(d.S1.BChAR.weightedC)

%% Inverse Gamma Parameters for each method

invGam_BMedian = [d.S1.BMedian.invGam.fitInfo.alpha, d.S1.BMedian.invGam.fitInfo.beta];
invGam_BChAR = [d.S1.BChAR.invGam.fitInfo.alpha, d.S1.BChAR.invGam.fitInfo.beta];
invGam_RSR0AR = [d.S1.RSR0AR.invGam.fitInfo.alpha, d.S1.RSR0AR.invGam.fitInfo.beta];
invGam_RSR500AR = [d.S1.RSR500AR.invGam.fitInfo.alpha, d.S1.RSR500AR.invGam.fitInfo.beta];
invGam_RSR1000AR =[d.S1.RSR1000AR.invGam.fitInfo.alpha, d.S1.RSR1000AR.invGam.fitInfo.beta];

invGam_parameters = [invGam_BMedian; invGam_BChAR; invGam_RSR0AR; invGam_RSR500AR; invGam_RSR1000AR]

%% Gamma parameters for each method

Gam_BMedian = [d.S1.BMedian.Gam.fitInfo.alpha, d.S1.BMedian.Gam.fitInfo.beta];
Gam_BChAR = [d.S1.BChAR.Gam.fitInfo.alpha, d.S1.BChAR.Gam.fitInfo.beta];
Gam_RSR0AR = [d.S1.RSR0AR.Gam.fitInfo.alpha, d.S1.RSR0AR.Gam.fitInfo.beta];
Gam_RSR500AR = [d.S1.RSR500AR.Gam.fitInfo.alpha, d.S1.RSR500AR.Gam.fitInfo.beta];
Gam_RSR1000AR =[d.S1.RSR1000AR.Gam.fitInfo.alpha, d.S1.RSR1000AR.Gam.fitInfo.beta];

Gam_parameters = [Gam_BMedian; Gam_BChAR; Gam_RSR0AR; Gam_RSR500AR; Gam_RSR1000AR]

%% LN parameters for each method

LN_BMedian = [d.S1.BMedian.LN.fitInfo.mu, d.S1.BMedian.LN.fitInfo.Sigma];
LN_BChAR = [d.S1.BChAR.LN.fitInfo.mu, d.S1.BChAR.LN.fitInfo.Sigma];
LN_RSR0AR = [d.S1.RSR0AR.LN.fitInfo.mu, d.S1.RSR0AR.LN.fitInfo.Sigma];
LN_RSR500AR = [d.S1.RSR500AR.LN.fitInfo.mu, d.S1.RSR500AR.LN.fitInfo.Sigma];
LN_RSR1000AR =[d.S1.RSR1000AR.LN.fitInfo.mu, d.S1.RSR1000AR.LN.fitInfo.Sigma];

LN_parameters = [LN_BMedian; LN_BChAR; LN_RSR0AR; LN_RSR500AR; LN_RSR1000AR]

%%
% Provide evidence that fitting to the pooled data
% of many MC runs provides results that are representative of fitting to
% the data of each MC run

%Log Bin Edges
lbe = [-4:0.1:4];

IRd = d.S1.RSR500IR;
ARd = d.S1.RSR500AR;

%Plot pooled histograms, Fit to each MC run, Fit to pooled data
figure;
subplot(2,4,1)
ARvsIRplots_wph(ARd.LN, IRd.LN, ARd, lbe)
text(0.02, 0.9, 'A', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

subplot(2,4,2)
ARvsIRplots_wph(ARd.MLN, IRd.MLN, ARd, lbe)
text(0.02, 0.9, 'B', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

subplot(2,4,3)
ARvsIRplots_wph(ARd.Gam, IRd.Gam, ARd, lbe)
text(0.02, 0.9, 'C', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

subplot(2,4,4)
ARvsIRplots_wph(ARd.invGam, IRd.invGam, ARd, lbe)
text(0.02, 0.9, 'D', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')


subplot(2,4,5)
hold on
for i = 1:400;
plot(IRd.LN.fits{i}.mu, sqrt(IRd.LN.fits{i}.Sigma), '.k')
end
plot(ARd.LN.fitInfo.mu, sqrt(ARd.LN.fitInfo.Sigma), 'or', 'MarkerFaceColor', 'r')
xlabel("Mu")
ylabel("Sigma")
text(0.02, 0.9, 'E', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

subplot(2,4,6)
hold on
for i = 1:400;
plot(IRd.MLN.nSR.mu(i), sqrt(IRd.MLN.nSR.var(i)), '.k')
end
plot(ARd.MLN.nSR.mu, sqrt(ARd.MLN.nSR.var), 'or', 'MarkerFaceColor', 'r')
xlabel("Mean")
ylabel("Standard Deviation")
text(0.02, 0.9, 'F', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

subplot(2,4,7)
hold on
for i = 1:400;
    plot(IRd.Gam.fits{i}.alpha, IRd.Gam.fits{i}.beta, '.k')
end
plot(ARd.Gam.fitInfo.alpha, ARd.Gam.fitInfo.beta, 'or', 'MarkerFaceColor', 'r')
xlabel("Alpha")
ylabel("Beta")
text(0.02, 0.9, 'G', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

subplot(2,4,8)
hold on
for i=1:400;
    plot(IRd.invGam.fits{i}.alpha, IRd.invGam.fits{i}.beta, '.k')
end
plot(ARd.invGam.fitInfo.alpha, ARd.invGam.fitInfo.beta, 'or', 'MarkerFaceColor', 'r')
xlabel("Alpha")
ylabel("Beta")
text(0.02, 0.9, 'H', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')



function ARvsIRplots_wph(ARd, IRd, AR, logBinEdges)
% Plot AR fit on top of IR fit, with pooled histogram underneath
AR_logHC = histcounts(log(AR.weightedC), logBinEdges);
normF = 1./sum((AR_logHC.*uniquetol(diff(logBinEdges), 1e-5)));
hold on
histogram('BinCounts', AR_logHC, 'BinEdges', logBinEdges, 'HandleVisibility', 'off', 'FaceColor', '[0.8 0.8 0.8]')
plot(IRd.lnSR.x, IRd.lnSR.px(:,1).*(1/normF), 'Color', [0 0 0 0.1], 'DisplayName', 'Individual Run MLN')
plot(IRd.lnSR.x, IRd.lnSR.px(:,2:end).*(1/normF), 'Color', [0 0 0 0.1], 'HandleVisibility', 'off')
plot(ARd.lnSR.x, ARd.lnSR.px.*(1/normF), 'Color', "r", 'LineWidth', 1, 'DisplayName', 'MLN')
xlim([-2.5 2.5])
ylim([0 1.5e5])
xlabel("log(NSR)")
ylabel("Pooled Counts")
end


function ARvsIRplots_muvar(ARd, IRd)
% Plot mu/var of IR fits and AR fit
hold on
plot(IRd.lnSR.mu, IRd.lnSR.var, '.', 'Color', [0 0 0 0.1])
plot(ARd.lnSR.mu, ARd.lnSR.var, 'o', 'Color', [1 0 0], 'MarkerFaceColor', [1 0 0])
xlabel("mean log(NSR)")
ylabel("var log(NSR)")
end