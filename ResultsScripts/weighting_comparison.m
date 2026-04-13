%% Compare distributions for different weighting approaches
%Load depth weighting results, store as dw
load("../Results/dataT_All1_RLGtrue_BchronJun2_3Apr26_fit3Apr26_depthweight.mat")
dd = d;

%Load age weighting results, store as da
load("../Results/dataT_All1_RLGtrue_BchronJun2_3Apr26_fit3Apr26_ageweight.mat")
da = d;

%Load no weighting results,store as dn
load("../Results/dataT_All1_RLGtrue_BchronJun2_3Apr26_fit3Apr26_noweight.mat")
dn = d;

%% Set up histogram settings
%Create binedges
maxbinedge   = 20;
bw           = 0.1; %bin width
binEdges     =  0:bw:maxbinedge; %nSR
invBinEdges  =  0:bw:maxbinedge; %invnSR
logBinEdges  = -5:bw:5;          %lognSR
logBinCenters = (logBinEdges(1:end-1)+logBinEdges(2:end)).*0.5;

%Set text size
set(groot, 'DefaultAxesFontSize', 12)
set(groot, 'DefaultTextFontSize', 12)

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

%The counts in each of these histograms is equivalent to the number of
%counts * total weighting * weighting inflator * number of runs.

BIC_table_RSR500 = array2table(...
    round([...
        dd.S1.RSR500AR.LN.fitInfo.BICtaeheefix, dd.S1.RSR500AR.MLN.fitInfo.BICtaeheefix, dd.S1.RSR500AR.Gam.fitInfo.BICtaeheefix, dd.S1.RSR500AR.invGam.fitInfo.BICtaeheefix;
        da.S1.RSR500AR.LN.fitInfo.BICtaeheefix, da.S1.RSR500AR.MLN.fitInfo.BICtaeheefix, da.S1.RSR500AR.Gam.fitInfo.BICtaeheefix, da.S1.RSR500AR.invGam.fitInfo.BICtaeheefix;
        dn.S1.RSR500AR.LN.fitInfo.BICtaeheefix, dn.S1.RSR500AR.MLN.fitInfo.BICtaeheefix, dn.S1.RSR500AR.Gam.fitInfo.BICtaeheefix, dn.S1.RSR500AR.invGam.fitInfo.BICtaeheefix;]), ...
    'VariableNames', {'LN', 'MLN', 'Gam', 'invGam'}, ...
    'RowNames',      {'Depth', 'Age', 'None'})

%% Plot BSamp distribution for all

figure;
subplot(3,1,1)
depthw_logHC = histcounts(log(dd.S1.BSampAR.weightedC), logBinEdges);
normF = 1./sum((depthw_logHC.*uniquetol(diff(logBinEdges), 1e-5)));
histogram(log(dd.S1.BSampAR.weightedC), 'BinEdges', logBinEdges, 'FaceColor', '[0.8 0.8 0.8]')
hold on
plot(dd.S1.BSampAR.LN.lnSR.x, dd.S1.BSampAR.LN.lnSR.px*(1/normF), 'LineWidth', 1, 'DisplayName','LN', 'Color', 'r')
plot(dd.S1.BSampAR.MLN.lnSR.x, dd.S1.BSampAR.MLN.lnSR.px*(1/normF), 'LineWidth', 1, 'DisplayName','MLN', 'Color', H(5))
plot(dd.S1.BSampAR.Gam.lnSR.x, dd.S1.BSampAR.Gam.lnSR.px*(1/normF),  'LineWidth', 1, 'DisplayName','Gam', 'Color', H(4) )
plot(dd.S1.BSampAR.invGam.lnSR.x, dd.S1.BSampAR.invGam.lnSR.px*(1/normF),  'LineWidth', 1, 'DisplayName','invGam', 'Color', 'b')
title("Depth Weighted")
ylabel("Counts")
xlim([-3 3])
subplot(3,1,2)
agew_logHC = histcounts(log(da.S1.BSampAR.weightedC), logBinEdges);
normF = 1./sum((agew_logHC.*uniquetol(diff(logBinEdges), 1e-5)));
histogram(log(da.S1.BSampAR.weightedC), 'BinEdges', logBinEdges, 'FaceColor', '[0.8 0.8 0.8]')
hold on
plot(da.S1.BSampAR.LN.lnSR.x, da.S1.BSampAR.LN.lnSR.px*(1/normF), 'LineWidth', 1, 'DisplayName','LN', 'Color', 'r')
plot(da.S1.BSampAR.MLN.lnSR.x, da.S1.BSampAR.MLN.lnSR.px*(1/normF), 'LineWidth', 1, 'DisplayName','MLN', 'Color', H(5))
plot(da.S1.BSampAR.Gam.lnSR.x, da.S1.BSampAR.Gam.lnSR.px*(1/normF),  'LineWidth', 1, 'DisplayName','Gam', 'Color', H(4) )
plot(da.S1.BSampAR.invGam.lnSR.x, da.S1.BSampAR.invGam.lnSR.px*(1/normF),  'LineWidth', 1, 'DisplayName','invGam', 'Color', 'b')

title("Age Weighted")
ylabel("Counts")
xlim([-3 3])
subplot(3,1,3)
now_logHC = histcounts(log(dn.S1.BSampAR.weightedC), logBinEdges);
normF = 1./sum((now_logHC.*uniquetol(diff(logBinEdges), 1e-5)));
histogram(log(dn.S1.BSampAR.weightedC), 'BinEdges', logBinEdges, 'FaceColor', '[0.8 0.8 0.8]')
hold on
plot(dn.S1.BSampAR.LN.lnSR.x, dn.S1.BSampAR.LN.lnSR.px*(1/normF), 'LineWidth', 1, 'DisplayName','LN', 'Color', 'r')
plot(dn.S1.BSampAR.MLN.lnSR.x, dn.S1.BSampAR.MLN.lnSR.px*(1/normF), 'LineWidth', 1, 'DisplayName','MLN', 'Color', H(5))
plot(dn.S1.BSampAR.Gam.lnSR.x, dn.S1.BSampAR.Gam.lnSR.px*(1/normF),  'LineWidth', 1, 'DisplayName','Gam', 'Color', H(4) )
plot(dn.S1.BSampAR.invGam.lnSR.x, dn.S1.BSampAR.invGam.lnSR.px*(1/normF),  'LineWidth', 1, 'DisplayName','invGam', 'Color', 'b')
title("Not Weighted")
ylabel("Counts")
xlim([-3 3])
xlabel("log(NSR)")

BIC_table_BSamp = array2table(...
    round([...
        dd.S1.BSampAR.LN.fitInfo.BICtaeheefix, dd.S1.BSampAR.MLN.fitInfo.BICtaeheefix, dd.S1.BSampAR.Gam.fitInfo.BICtaeheefix, dd.S1.BSampAR.invGam.fitInfo.BICtaeheefix;
        da.S1.BSampAR.LN.fitInfo.BICtaeheefix, da.S1.BSampAR.MLN.fitInfo.BICtaeheefix, da.S1.BSampAR.Gam.fitInfo.BICtaeheefix, da.S1.BSampAR.invGam.fitInfo.BICtaeheefix;
        dn.S1.BSampAR.LN.fitInfo.BICtaeheefix, dn.S1.BSampAR.MLN.fitInfo.BICtaeheefix, dn.S1.BSampAR.Gam.fitInfo.BICtaeheefix, dn.S1.BSampAR.invGam.fitInfo.BICtaeheefix;]), ...
    'VariableNames', {'LN', 'MLN', 'Gam', 'invGam'}, ...
    'RowNames',      {'Depth', 'Age', 'None'})