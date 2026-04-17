%% Header
%This script produces the figures and tables that are relevant for the
%manuscript.

%Figure 1 - Map of radiocarbon dated cores used, as well as histogram of
%seafloor depths at core sites and average SR across the core.

%Table 1 - Total amount of sediment and number of age pairs included in 
% all cores analyzed for each sampling method

%Figure 2 - Histograms of weighted NSR (log scale) and \Delta t for 
%the BIGMACS results and the BMedian from this study.

%Figure 3 - Weighted histograms of NSR and Δt for BSamp (A, B),
%  RSR0 (C, D), RSR500 (E, F) and RSR1000 (G, H). To represent all 
% iterations, the gray bars show the median value of each bin and the 
% error bars show bounds of the 68% confidence intervals from each bin

%Table 2 - BIC of the fitted distributions for each sampling method used.

%Figure 4 - Panel A shows the weighted histogram of pooled data from 
% RSR500 method, with the 4 fitted distributions 
% (Lognormal – LN; Mixture Lognormal – MLN; Gamma – G; Inverse Gamma – IG). 
% Panel B shows the same weighted histogram of pooled data from RSR500 
% run with the lognormal fit to the pooled data (red) and the lognormal 
% fits of each individual RSR500 run (thin gray lines).

% Supplementary Figure A1 - RSR500 pooled histogram, lognormal and inverse
% gamma, plotted on inverse SR axis.

% Supplementary Figure A2 - 

%% Add necessary paths
addpath("../Functions/")

%Set text size for figures
set(groot, 'DefaultAxesFontSize', 12)
set(groot, 'DefaultTextFontSize', 12)

%% Set up some general histogram settings
%Create binedges
maxbinedge   = 20;
bw           = 0.1; %bin width
binEdges     =  0:bw:maxbinedge; %nSR
invBinEdges  =  0:bw:maxbinedge; %invnSR
logBinEdges  = -5:bw:5;          %lognSR
logBinCenters = (logBinEdges(1:end-1)+logBinEdges(2:end)).*0.5;

%% Input and Analyse BIGMACS data
%Load BIGMACS files
BMlognorm = readtable("../BIGMACSdata/lognormal.txt");
BM.nSR.x = BMlognorm.Var1';
BM.nSR.px = BMlognorm.Var2;
BM.weightedC = readmatrix("../BIGMACSdata/Lin2014_sedrateratio_cm_wo_NaN.txt"); % all weighted counts
BM.weightedC = BM.weightedC(:,4);
BM.sedLength = length(BM.weightedC); % in cm; Because every NSR value represents 1cm of sediment
BM.sedTimeSpan = 544; % kyr; From Lin et al. 2014;
BM.numCpairs = 343; % known from running code used for Lin2014 (not available to others)

BM.TM = readmatrix("../BIGMACSdata/transition_parameter.txt");
[MLN_BIGMACS, ~, BM.gm2fit] = fitMixLogNorm(BM.weightedC, BM.nSR.x, 2, 3, 343);
[LN_BIGMACS, ~, BM.gm1fit] = fitMixLogNorm(BM.weightedC, BM.nSR.x, 1, 3, 343);
[~, invGam_BIGMACS, BM.invgamfit] = fitInvGamma(BM.weightedC, BM.nSR.x, 343);
[~,BM.gamfit] = fitGamma(BM.weightedC, BM.nSR.x, 343);
[BM.lnSR.x, BM.lnSR.px] = px_to_pfx(BM.nSR.x, BM.nSR.px, @log);
BMpdf.numParams = 5;
BMpdf.x = BM.lnSR.x;
BMpdf.px = BM.lnSR.px;

%% Load data and fits
%load("../Results/dataT_All1_RLGtrue_BchronJun2_3Apr26_fit3Apr26_depthweight.mat");
load("../Results/dataT_All1_RLGtrue_Bchron8Apr26_fit9Apr26_depthweight");

%% Take data of interest
%----- Find corenames, lats, longs, depths of cores included
dataTbl = d.dataT(d.S1.chooseLog, :);
%----- Make map with locations denoted as red squares

% Define tick locations
lonTicks = -180:20:180;
latTicks = -40:20:40;
tickLen = 2; % in degrees

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

[lonTick_x, lonTick_y] = projfwd(proj, repmat(latTicks(1), size(lonTicks)), lonTicks);
[latTick_x, latTick_y] = projfwd(proj, latTicks, repmat(lonTicks(1), size(latTicks)));

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
tightmap

%Manipulate y tick labels
yticklabelHandles = findall(gcf,'Type','text','Tag','PLabel');
nyLabeles = numel(yticklabelHandles);
ytickPosition = cell2mat(get(yticklabelHandles,'Position'));
ytickPosition(:,1) = ytickPosition(:,1) + range(ax.XLim)*.02;
set(yticklabelHandles, {'Position'}, mat2cell(ytickPosition, ones(nyLabeles,1),3))
tightmap

% Histogram of core depths
subplot(2,2,3)
histogram(dataTbl.depths./1000, 0:0.5:6, 'FaceColor', [0.8 0.8 0.8])
xlim([0 6])
xlabel("Depth (km)")
ylabel("Counts")
xticks(0:1:6)
ylims = ylim; ylim(ylims.*1.1)
set(gca, 'Position', [0.13 0.25 0.33 0.25])
text(0.02, 0.9, 'B', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

% Histogram of core average SRs
subplot(2,2,4)
histogram(dataTbl.meanSR, 0:5:90, 'FaceColor', [0.8 0.8 0.8])
xlabel('Average SR (cm/kyr)')
ylabel('Counts')
xlim([0 90])
ylims = ylim; ylim(ylims.*1.1)
set(gca, 'Position', [0.57 0.25 0.33 0.25])
text(0.02, 0.9, 'C', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')


%% %Create table of useful information
dStrus = {BM, d.S1.BMedian, d.S1.BSampAR, d.S1.RSR0AR, d.S1.RSR500AR, d.S1.RSR1000AR};
dStrusStrings = ["BIGMACS", "BMedian", "BSamp", "RSR0", "RSR500", "RSR1000"];

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
        MeanSedLength(i,1) = dStru.sedLength;
        MeanSedTime(i,1)   = dStru.sedTimeSpan;
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

%% Table 1 : Amount of data/sediment used in each method
dA.summT(:, 1:3) 

%% Figure 2: Plot of BIGMACS and BMedian NSR with normalized BIGMACS pdfs and histograms of dt

%Load NSR information of BIGMACS data
load("../Results/Lin2014_dts.mat")

%Get weighted histogram of dts for BIGMACS
BM_dt_weightedBC = makeWeightedBinCounts(asr(:,2), asr(:,3), 0:0.1:10);

%Set up figure
fA = figure;
fA.Position = [100, 100, 1000, 400]; % Wide and short: [left, bottom, width, height]

%Plot histogram of BIGMACS weighted NSR counts
subplot(2,2,1)
BM.logHC = histcounts(log(BM.weightedC), logBinEdges);
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

%Plot histogram of BIGMACS weighted dts
subplot(2,2,2)
hold on
histogram('BinCounts', BM_dt_weightedBC, 'BinEdges', 0:0.1:10,  'FaceColor', [0.8 0.8 0.8])
text(0.02, 0.9, 'B', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')
box on
ylabel("cm")
xlabel("\Deltat (kyr)")
xlim([0 5])
ylim([0 2000])

%Plot histogram of BMedian weighted NSR counts
subplot(2,2,3)
BMedian_logHC = histcounts(log(d.S1.BMedian.weightedC), logBinEdges)./d.S1.fitS.PooledRuns.weightRepInflator;
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

%Plot histogram of BMedian weighted dts
subplot(2, 2, 4)
histogram('BinCounts', d.S1.BMedian.agediffsWC, 'BinEdges', 0:0.1:10,  'FaceColor', [0.8 0.8 0.8])
text(0.02, 0.9, 'D', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')
box on
ylabel("cm")
xlabel("\Deltat (kyr)")
xlim([0 5])
ylim([0 2000])


%% Figure 3 All sampling approaches histogram of median with 68th percentile bars + dt histogram


numruns = length(d.S1.BSampIR.OneRunDatas);
fB = figure;
fB.Position = [100, 100, 1000, 800]; % Wide and short: [left, bottom, width, height]
subplot(4,2,1)
A = sort(d.S1.BSampIR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', A(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceColor', '[0.8 0.8 0.8]')
hold on
%plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds', 'Color', 'k')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
xlabel("log(NSR)")
ylabel("cm")
text(0.02, 0.9, 'A', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

subplot(4,2,2)
B = sort(d.S1.BSampIR.agediffsWbC', 1);
hold on
box on
histogram('BinCounts', B(numruns*0.5, :), 'BinEdges', 0:0.1:10, 'FaceColor', '[0.8 0.8 0.8]')
errorbar(0.05:0.1:9.95, B(numruns*0.5, :),(B(numruns*0.5, :)-B(numruns*0.16, :)),(B(numruns*0.84, :)-B(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds', 'Color', 'k')
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
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds', 'Color', 'k')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])

ylabel("cm")
xlabel("log(NSR)")
text(0.02, 0.9, 'C', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

subplot(4,2,4)
B = sort(d.S1.RSR0IR.agediffsWbC', 1);
hold on
box on
histogram('BinCounts', B(numruns*0.5, :), 'BinEdges', 0:0.1:10, 'FaceColor', '[0.8 0.8 0.8]')
errorbar(0.05:0.1:9.95, B(numruns*0.5, :),(B(numruns*0.5, :)-B(numruns*0.16, :)),(B(numruns*0.84, :)-B(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds', 'Color', 'k')
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
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds', 'Color', 'k')
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
errorbar(0.05:0.1:9.95, B(numruns*0.5, :),(B(numruns*0.5, :)-B(numruns*0.16, :)),(B(numruns*0.84, :)-B(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds', 'Color', 'k')
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
errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds', 'Color', 'k')
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
errorbar(0.05:0.1:9.95, B(numruns*0.5, :),(B(numruns*0.5, :)-B(numruns*0.16, :)),(B(numruns*0.84, :)-B(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds', 'Color', 'k')
xlim([0 6])
ylim([0 2000])
xlabel("\Deltat (kyr)")
ylabel("cm")
text(0.02, 0.9, 'H', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

%Make all the y limits the same for first column of plots
max_ylim = 0;
for i = [1,3,5,7]
    subplot(4,2,i)
    new_ylim = max(ylim);
    if new_ylim > max_ylim
        max_ylim = new_ylim;
    end
end

for i = [1,3,5,7]
    subplot(4,2,i)
    ylim([0 new_ylim.*1.05])
end

%Make all the y limits the same for second column of plots
max_ylim = 0;
for i = [2,4,6,8]
    subplot(4,2,i)
    new_ylim = max(ylim);
    if new_ylim > max_ylim
        max_ylim = new_ylim;
    end
end

for i = [2,4,6,8]
    subplot(4,2,i)
    ylim([0 new_ylim.*1.05])
end

%% Create Table 2 showing BIC results

%This is the BIC (with the correction to the likelihood) for each 

BIC_table = array2table(...
    round([...
        d.S1.BMedian.LN.fitInfo.BICtaeheefix,   d.S1.BMedian.MLN.fitInfo.BICtaeheefix,   d.S1.BMedian.Gam.fitInfo.BICtaeheefix,   d.S1.BMedian.invGam.fitInfo.BICtaeheefix;
        d.S1.BSampAR.LN.fitInfo.BICtaeheefix,   d.S1.BSampAR.MLN.fitInfo.BICtaeheefix,   d.S1.BSampAR.Gam.fitInfo.BICtaeheefix,   d.S1.BSampAR.invGam.fitInfo.BICtaeheefix;
        d.S1.RSR0AR.LN.fitInfo.BICtaeheefix,    d.S1.RSR0AR.MLN.fitInfo.BICtaeheefix,    d.S1.RSR0AR.Gam.fitInfo.BICtaeheefix,    d.S1.RSR0AR.invGam.fitInfo.BICtaeheefix;
        d.S1.RSR500AR.LN.fitInfo.BICtaeheefix,  d.S1.RSR500AR.MLN.fitInfo.BICtaeheefix,  d.S1.RSR500AR.Gam.fitInfo.BICtaeheefix,  d.S1.RSR500AR.invGam.fitInfo.BICtaeheefix;
        d.S1.RSR1000AR.LN.fitInfo.BICtaeheefix, d.S1.RSR1000AR.MLN.fitInfo.BICtaeheefix, d.S1.RSR1000AR.Gam.fitInfo.BICtaeheefix, d.S1.RSR1000AR.invGam.fitInfo.BICtaeheefix]), ...
    'VariableNames', {'LN', 'MLN', 'Gam', 'invGam'}, ...
    'RowNames',      {'BMedian', 'BSamp', 'RSR0', 'RSR500', 'RSR1000'})


%% Figure 4 - Plot of pooled samples + fits

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
ylims = ylim;
text(0.02, 0.9, 'A', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')


%And the individual runs as well
subplot(2,1,2)
hold on;
AR_logHC = histcounts(log(AR.weightedC), logBinEdges);
normF = 1./sum((AR_logHC.*uniquetol(diff(logBinEdges), 1e-5)));
histogram('BinCounts', AR_logHC, 'BinEdges', logBinEdges, 'HandleVisibility', 'off', 'FaceColor', '[0.8 0.8 0.8]')
plot(IR.LN.lnSR.x, IR.LN.lnSR.px(:,1).*(1/normF), 'Color', [0 0 0 0.2], 'DisplayName', 'Individual Run LN')
plot(IR.LN.lnSR.x, IR.LN.lnSR.px(:,2:end).*(1/normF), 'Color', [0 0 0 0.05], 'HandleVisibility', 'off')
plot(AR.LN.lnSR.x, AR.LN.lnSR.px.*(1/normF), 'Color', "r", 'LineWidth', 1, 'DisplayName', 'LN')
%plot(BM.lnSR.x, BM.lnSR.px.*(1/normF), 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'DisplayName', 'BIGMACS')
xlim([-2.5 2.5])
ylim(ylims)
xlabel("log(NSR)")
ylabel("Pooled Counts")
legend()
text(0.02, 0.9, 'B', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')


%% Supplementary Plots and Tables
%% Plot of pooled samples + fits on inverse NSR axis

%Choose sampling method of interest
AR = d.S1.RSR500AR;

%define histogram bins to use
invBinEdges = 0:0.1:6;

%Create pdfs on inverse SR
f_inv = @(x)(1./x);
[AR.LN.invSR.x, AR.LN.invSR.px] = px_to_pfx(AR.LN.nSR.x, AR.LN.nSR.px, f_inv);
[AR.invGam.invSR.x, AR.invGam.invSR.px] = px_to_pfx(AR.invGam.nSR.x, AR.invGam.nSR.px, f_inv);

%Create the bacon default distribution
shape = 1.5;
bacon_invGam = gampdf(AR.LN.invSR.x, shape, 1/shape);

figure;
hold on
AR_invHC = histcounts(1./AR.weightedC, invBinEdges);
normF = 1./sum((AR_invHC.*uniquetol(diff(invBinEdges), 1e-5)));
histogram('BinCounts', AR_invHC, 'BinEdges', invBinEdges, 'HandleVisibility', 'off', 'FaceColor', '[0.8 0.8 0.8]')
plot(AR.LN.invSR.x, AR.LN.invSR.px*(1/normF), 'r', 'LineWidth', 1, 'DisplayName','LN')
plot(AR.invGam.invSR.x, AR.invGam.invSR.px*(1/normF),  'LineWidth', 1, 'DisplayName','IG', 'Color', 'b')
plot(AR.LN.invSR.x, bacon_invGam*(1/normF),  'LineWidth', 1, 'DisplayName','Bacon Default', 'Color', 'g')
xlim([0 5])
xlabel("inverse NSR")
ylabel("Pooled Counts")
legend()
%ylim([0 1.2e5])
ylims = ylim;


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

%Give all similar figures the same ylims
max_ylim = 0;
for i = [1,2,3,4]
    subplot(2,4,i)
    new_ylim = max(ylim);
    if new_ylim > max_ylim
        max_ylim = new_ylim;
    end
end

for i = [1,2,3,4]
    subplot(2,4,i)
    ylim([0 new_ylim.*1.05])
end




%% Compare distributions for different weighting approaches
%Load age weighting results, store as da
load("../Results/dataT_All1_RLGtrue_BchronJun2_3Apr26_fit3Apr26_ageweight.mat")
da = d;

%Load no weighting results,store as dn
load("../Results/dataT_All1_RLGtrue_BchronJun2_3Apr26_fit3Apr26_noweight.mat")
dn = d;

%Load depth weighting results, store as dw
load("../Results/dataT_All1_RLGtrue_BchronJun2_3Apr26_fit3Apr26_depthweight.mat")
dd = d;
%% Plot RSR500 AR distribution for all

RGB = orderedcolors("gem");
H = rgb2hex(RGB);

figure;
subplot(3,1,1)
depthw_logHC = histcounts(log(dd.S1.RSR500AR.weightedC), logBinEdges);
normF = 1./sum((depthw_logHC.*uniquetol(diff(logBinEdges), 1e-5)));
histogram(log(dd.S1.RSR500AR.weightedC), 'BinEdges', logBinEdges, 'FaceColor', '[0.8 0.8 0.8]', 'HandleVisibility','off')
hold on
plot(dd.S1.RSR500AR.LN.lnSR.x, dd.S1.RSR500AR.LN.lnSR.px*(1/normF), 'LineWidth', 1, 'DisplayName','LN', 'Color', 'r')
plot(dd.S1.RSR500AR.MLN.lnSR.x, dd.S1.RSR500AR.MLN.lnSR.px*(1/normF), 'LineWidth', 1, 'DisplayName','MLN', 'Color', H(5))
plot(dd.S1.RSR500AR.Gam.lnSR.x, dd.S1.RSR500AR.Gam.lnSR.px*(1/normF),  'LineWidth', 1, 'DisplayName','G', 'Color', H(4) )
plot(dd.S1.RSR500AR.invGam.lnSR.x, dd.S1.RSR500AR.invGam.lnSR.px*(1/normF),  'LineWidth', 1, 'DisplayName','IG', 'Color', 'b')
legend()
%title("Depth Weighted")
text(0.02, 0.85, 'A', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')
ylabel("Counts")
xlim([-3 3])
subplot(3,1,2)
agew_logHC = histcounts(log(da.S1.RSR500AR.weightedC), logBinEdges);
normF = 1./sum((agew_logHC.*uniquetol(diff(logBinEdges), 1e-5)));
histogram(log(da.S1.RSR500AR.weightedC), 'BinEdges', logBinEdges, 'FaceColor', '[0.8 0.8 0.8]')
hold on
plot(da.S1.RSR500AR.LN.lnSR.x, da.S1.RSR500AR.LN.lnSR.px*(1/normF), 'LineWidth', 1, 'DisplayName','LN', 'Color', 'r')
plot(da.S1.RSR500AR.MLN.lnSR.x, da.S1.RSR500AR.MLN.lnSR.px*(1/normF), 'LineWidth', 1, 'DisplayName','MLN', 'Color', H(5))
plot(da.S1.RSR500AR.Gam.lnSR.x, da.S1.RSR500AR.Gam.lnSR.px*(1/normF),  'LineWidth', 1, 'DisplayName','G', 'Color', H(4) )
plot(da.S1.RSR500AR.invGam.lnSR.x, da.S1.RSR500AR.invGam.lnSR.px*(1/normF),  'LineWidth', 1, 'DisplayName','IG', 'Color', 'b')
text(0.02, 0.85, 'B', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')
%title("Age Weighted")
ylabel("Counts")
xlim([-3 3])
subplot(3,1,3)
now_logHC = histcounts(log(dn.S1.RSR500AR.weightedC), logBinEdges);
normF = 1./sum((now_logHC.*uniquetol(diff(logBinEdges), 1e-5)));
histogram(log(dn.S1.RSR500AR.weightedC), 'BinEdges', logBinEdges, 'FaceColor', '[0.8 0.8 0.8]')
hold on
plot(dn.S1.RSR500AR.LN.lnSR.x, dn.S1.RSR500AR.LN.lnSR.px*(1/normF), 'LineWidth', 1, 'DisplayName','LN', 'Color', 'r')
plot(dn.S1.RSR500AR.MLN.lnSR.x, dn.S1.RSR500AR.MLN.lnSR.px*(1/normF), 'LineWidth', 1, 'DisplayName','MLN', 'Color', H(5))
plot(dn.S1.RSR500AR.Gam.lnSR.x, dn.S1.RSR500AR.Gam.lnSR.px*(1/normF),  'LineWidth', 1, 'DisplayName','G', 'Color', H(4) )
plot(dn.S1.RSR500AR.invGam.lnSR.x, dn.S1.RSR500AR.invGam.lnSR.px*(1/normF),  'LineWidth', 1, 'DisplayName','IG', 'Color', 'b')
%title("Not Weighted")
text(0.02, 0.85, 'C', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')
ylabel("Counts")
xlim([-3 3])
xlabel("log(NSR)")

%Get BIC table of each
BIC_table_RSR500 = array2table(...
    round([...
        dd.S1.RSR500AR.LN.fitInfo.BICtaeheefix, dd.S1.RSR500AR.MLN.fitInfo.BICtaeheefix, dd.S1.RSR500AR.Gam.fitInfo.BICtaeheefix, dd.S1.RSR500AR.invGam.fitInfo.BICtaeheefix;
        da.S1.RSR500AR.LN.fitInfo.BICtaeheefix, da.S1.RSR500AR.MLN.fitInfo.BICtaeheefix, da.S1.RSR500AR.Gam.fitInfo.BICtaeheefix, da.S1.RSR500AR.invGam.fitInfo.BICtaeheefix;
        dn.S1.RSR500AR.LN.fitInfo.BICtaeheefix, dn.S1.RSR500AR.MLN.fitInfo.BICtaeheefix, dn.S1.RSR500AR.Gam.fitInfo.BICtaeheefix, dn.S1.RSR500AR.invGam.fitInfo.BICtaeheefix;]), ...
    'VariableNames', {'LN', 'MLN', 'Gam', 'invGam'}, ...
    'RowNames',      {'Depth', 'Age', 'None'})

%% Check influence of not normalizing data
load("../Results/dataT_All1_RLGtrue_DS0p05_Dec9_fit26Mar26_noweight_SR2_400R")
dSR = d;

BIC_table2 = array2table(...
    round([...
        d.S1.BMedian.LN.fitInfo.BICtaeheefix,   d.S1.BMedian.MLN.fitInfo.BICtaeheefix,   d.S1.BMedian.Gam.fitInfo.BICtaeheefix,   d.S1.BMedian.invGam.fitInfo.BICtaeheefix;
        d.S1.BChAR.LN.fitInfo.BICtaeheefix,   d.S1.BChAR.MLN.fitInfo.BICtaeheefix,   d.S1.BChAR.Gam.fitInfo.BICtaeheefix,   d.S1.BChAR.invGam.fitInfo.BICtaeheefix;
        d.S1.New0AR.LN.fitInfo.BICtaeheefix,    d.S1.New0AR.MLN.fitInfo.BICtaeheefix,    d.S1.New0AR.Gam.fitInfo.BICtaeheefix,    d.S1.New0AR.invGam.fitInfo.BICtaeheefix;
        d.S1.New500AR.LN.fitInfo.BICtaeheefix,  d.S1.New500AR.MLN.fitInfo.BICtaeheefix,  d.S1.New500AR.Gam.fitInfo.BICtaeheefix,  d.S1.New500AR.invGam.fitInfo.BICtaeheefix;
        d.S1.New1000AR.LN.fitInfo.BICtaeheefix, d.S1.New1000AR.MLN.fitInfo.BICtaeheefix, d.S1.New1000AR.Gam.fitInfo.BICtaeheefix, d.S1.New1000AR.invGam.fitInfo.BICtaeheefix]), ...
    'VariableNames', {'LN', 'MLN', 'Gam', 'invGam'}, ...
    'RowNames',      {'BMedian', 'BCh', 'New0', 'New500', 'New1000'})

%% Helper functions

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