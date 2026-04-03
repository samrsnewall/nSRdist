%This script looks at the BIC of fits for model selection

% Add necessary paths
%addpath("../Functions/")

% Load Dataa
%dA = load("../Results/dataT_All1_RLGtrue_DS0p05_Dec9_fitDec9_depthweight_400R.mat");
binEdges = 0:0.1:10;

% distNames = ["LN", "MLN", "Gamma", "InvGamma"];
% resultNames = ["BMedian", "BSamp", "RSR0", "RSR500", "RSR1000"];
% BICTable = table('VariableNames',distNames, 'RowNames', resultNames);
 BMedianBICS = [d.S1.BMedian.LN.fits{1}.BICtaeheefix, d.S1.BMedian.MLN.fits{1}.BICtaeheefix, d.S1.BMedian.Gam.fits{1}.BICtaeheefix, d.S1.BMedian.invGam.fits{1}.BICtaeheefix];
% LNBICS = [d.S1.BMedian.LN.fits{1}.BICtaeheefix; d.S1.BMedian.MLN.fits{1}.BICtaeheefix; d.S1.BMedian.Gam.fits{1}.BICtaeheefix; d.S1.BMedian.invGam.fits{1}.BICtaeheefix];

%BMedian Fits
figure; 
subplot(2,2,1);
yyaxis left; histogram(d.S1.BMedian.weightedC{1}, 'BinEdges', binEdges);
yyaxis right; plot(d.S1.BMedian.MLN.nSR.x, d.S1.BMedian.MLN.nSR.px, 'LineWidth', 1);
title("BIC = "+num2str(d.S1.BMedian.MLN.fits{1}.BICtaeheefix, 3));
xlim([0 10]); ylim([0 1])
subplot(2,2,2);
yyaxis left; histogram(d.S1.BMedian.weightedC{1}, 'BinEdges', binEdges);
yyaxis right; plot(d.S1.BMedian.LN.nSR.x, d.S1.BMedian.LN.nSR.px, 'LineWidth', 1);
title("LN: BIC = "+num2str(d.S1.BMedian.LN.fits{1}.BICtaeheefix, 3));
xlim([0 10]); ylim([0 1])
subplot(2,2,3);
yyaxis left; histogram(d.S1.BMedian.weightedC{1}, 'BinEdges', binEdges);
yyaxis right; plot(d.S1.BMedian.Gam.nSR.x, d.S1.BMedian.Gam.nSR.px, 'LineWidth', 1);
title("Gam: BIC = "+num2str(d.S1.BMedian.Gam.fits{1}.BICtaeheefix,3));
xlim([0 10]); ylim([0 1])
subplot(2,2,4);
yyaxis left; histogram(d.S1.BMedian.weightedC{1}, 'BinEdges', binEdges);
yyaxis right; plot(d.S1.BMedian.invGam.nSR.x, d.S1.BMedian.invGam.nSR.px, 'LineWidth', 1);
title("invGam: BIC = "+num2str(d.S1.BMedian.invGam.fits{1}.BICtaeheefix,3));
xlim([0 10]); ylim([0 1])

%BSampAR fits
figure; 
subplot(2,2,1);
yyaxis left; histogram(d.S1.BSampAR.weightedC, 'BinEdges', binEdges);
yyaxis right; plot(d.S1.BSampAR.MLN.nSR.x, d.S1.BSampAR.MLN.nSR.px, 'LineWidth', 1);
title("BIC = "+num2str(d.S1.BSampAR.MLN.fitInfo.BICtaeheefix, 3));
xlim([0 10]); ylim([0 1])
subplot(2,2,2);
yyaxis left; histogram(d.S1.BSampAR.weightedC, 'BinEdges', binEdges);
yyaxis right; plot(d.S1.BSampAR.LN.nSR.x, d.S1.BSampAR.LN.nSR.px, 'LineWidth', 1);
title("LN: BIC = "+num2str(d.S1.BSampAR.LN.fitInfo.BICtaeheefix, 3));
xlim([0 10]); ylim([0 1])
subplot(2,2,3);
yyaxis left; histogram(d.S1.BSampAR.weightedC, 'BinEdges', binEdges);
yyaxis right; plot(d.S1.BSampAR.Gam.nSR.x, d.S1.BSampAR.Gam.nSR.px, 'LineWidth', 1);
title("Gam: BIC = "+num2str(d.S1.BSampAR.Gam.fitInfo.BICtaeheefix,3));
xlim([0 10]); ylim([0 1])
subplot(2,2,4);
yyaxis left; histogram(d.S1.BSampAR.weightedC, 'BinEdges', binEdges);
yyaxis right; plot(d.S1.BSampAR.invGam.nSR.x, d.S1.BSampAR.invGam.nSR.px, 'LineWidth', 1);
title("invGam: BIC = "+num2str(d.S1.BSampAR.invGam.fitInfo.BICtaeheefix,3));
xlim([0 10]); ylim([0 1])

%RSR0 fits
figure; 
subplot(2,2,1);
yyaxis left; histogram(log(d.S1.New0AR.weightedC));
yyaxis right; plot(d.S1.New0AR.MLN.lnSR.x, d.S1.New0AR.MLN.lnSR.px, 'LineWidth', 1);
title("BIC = "+num2str(d.S1.New0AR.MLN.fitInfo.BICtaeheefix, 3));
%xlim([0 10]); ylim([0 1])
subplot(2,2,2);
yyaxis left; histogram(log(d.S1.New0AR.weightedC));
yyaxis right; plot(d.S1.New0AR.LN.lnSR.x, d.S1.New0AR.LN.lnSR.px, 'LineWidth', 1);
title("LN: BIC = "+num2str(d.S1.New0AR.LN.fitInfo.BICtaeheefix, 3));
%xlim([0 10]); ylim([0 1])
subplot(2,2,3);
yyaxis left; histogram(log(d.S1.New0AR.weightedC));
yyaxis right; plot(d.S1.New0AR.Gam.lnSR.x, d.S1.New0AR.Gam.lnSR.px, 'LineWidth', 1);
title("Gam: BIC = "+num2str(d.S1.New0AR.Gam.fitInfo.BICtaeheefix,3));
%xlim([0 10]); ylim([0 1])
subplot(2,2,4);
yyaxis left; histogram(log(d.S1.New0AR.weightedC));
yyaxis right; plot(d.S1.New0AR.invGam.lnSR.x, d.S1.New0AR.invGam.lnSR.px, 'LineWidth', 1);
title("invGam: BIC = "+num2str(d.S1.New0AR.invGam.fitInfo.BICtaeheefix,3));
%xlim([0 10]); ylim([0 1])

%RSR500 fits
figure; 
subplot(2,2,1);
yyaxis left; histogram(log(d.S1.New500AR.weightedC));
yyaxis right; plot(d.S1.New500AR.MLN.lnSR.x, d.S1.New500AR.MLN.lnSR.px, 'LineWidth', 1);
title("BIC = "+num2str(d.S1.New500AR.MLN.fitInfo.BICtaeheefix, 3));
%xlim([0 10]); ylim([0 1])
subplot(2,2,2);
yyaxis left; histogram(log(d.S1.New500AR.weightedC));
yyaxis right; plot(d.S1.New500AR.LN.lnSR.x, d.S1.New500AR.LN.lnSR.px, 'LineWidth', 1);
title("LN: BIC = "+num2str(d.S1.New500AR.LN.fitInfo.BICtaeheefix, 3));
%xlim([0 10]); ylim([0 1])
subplot(2,2,3);
yyaxis left; histogram(log(d.S1.New500AR.weightedC));
yyaxis right; plot(d.S1.New500AR.Gam.lnSR.x, d.S1.New500AR.Gam.lnSR.px, 'LineWidth', 1);
title("Gam: BIC = "+num2str(d.S1.New500AR.Gam.fitInfo.BICtaeheefix,3));
%xlim([0 10]); ylim([0 1])
subplot(2,2,4);
yyaxis left; histogram(log(d.S1.New500AR.weightedC));
yyaxis right; plot(d.S1.New500AR.invGam.lnSR.x, d.S1.New500AR.invGam.lnSR.px, 'LineWidth', 1);
title("invGam: BIC = "+num2str(d.S1.New500AR.invGam.fitInfo.BICtaeheefix,3));
%xlim([0 10]); ylim([0 1])

%RSR1000 fits
figure; 
subplot(2,2,1);
yyaxis left; histogram(log(d.S1.New1000AR.weightedC));
yyaxis right; plot(d.S1.New1000AR.MLN.lnSR.x, d.S1.New1000AR.MLN.lnSR.px, 'LineWidth', 1);
title("BIC = "+num2str(d.S1.New1000AR.MLN.fitInfo.BICtaeheefix, 3));
%xlim([0 10]); ylim([0 1])
subplot(2,2,2);
yyaxis left; histogram(log(d.S1.New1000AR.weightedC));
yyaxis right; plot(d.S1.New1000AR.LN.lnSR.x, d.S1.New1000AR.LN.lnSR.px, 'LineWidth', 1);
title("LN: BIC = "+num2str(d.S1.New1000AR.LN.fitInfo.BICtaeheefix, 3));
%xlim([0 10]); ylim([0 1])
subplot(2,2,3);
yyaxis left; histogram(log(d.S1.New1000AR.weightedC));
yyaxis right; plot(d.S1.New1000AR.Gam.lnSR.x, d.S1.New1000AR.Gam.lnSR.px, 'LineWidth', 1);
title("Gam: BIC = "+num2str(d.S1.New1000AR.Gam.fitInfo.BICtaeheefix,3));
%xlim([0 10]); ylim([0 1])
subplot(2,2,4);
yyaxis left; histogram(log(d.S1.New1000AR.weightedC));
yyaxis right; plot(d.S1.New1000AR.invGam.lnSR.x, d.S1.New1000AR.invGam.lnSR.px, 'LineWidth', 1);
title("invGam: BIC = "+num2str(d.S1.New1000AR.invGam.fitInfo.BICtaeheefix,3));
%xlim([0 10]); ylim([0 1])

%RSR1500 fits
figure; 
subplot(2,2,1);
yyaxis left; histogram(log(d.S1.New1500AR.weightedC));
yyaxis right; plot(d.S1.New1500AR.MLN.lnSR.x, d.S1.New1500AR.MLN.lnSR.px, 'LineWidth', 1);
title("BIC = "+num2str(d.S1.New1500AR.MLN.fitInfo.BICtaeheefix, 3));
%xlim([0 10]); ylim([0 1])
subplot(2,2,2);
yyaxis left; histogram(log(d.S1.New1500AR.weightedC));
yyaxis right; plot(d.S1.New1500AR.LN.lnSR.x, d.S1.New1500AR.LN.lnSR.px, 'LineWidth', 1);
title("LN: BIC = "+num2str(d.S1.New1500AR.LN.fitInfo.BICtaeheefix, 3));
%xlim([0 10]); ylim([0 1])
subplot(2,2,3);
yyaxis left; histogram(log(d.S1.New1500AR.weightedC));
yyaxis right; plot(d.S1.New1500AR.Gam.lnSR.x, d.S1.New1500AR.Gam.lnSR.px, 'LineWidth', 1);
title("Gam: BIC = "+num2str(d.S1.New1500AR.Gam.fitInfo.BICtaeheefix,3));
%xlim([0 10]); ylim([0 1])
subplot(2,2,4);
yyaxis left; histogram(log(d.S1.New1500AR.weightedC));
yyaxis right; plot(d.S1.New1500AR.invGam.lnSR.x, d.S1.New1500AR.invGam.lnSR.px, 'LineWidth', 1);
title("invGam: BIC = "+num2str(d.S1.New1500AR.invGam.fitInfo.BICtaeheefix,3));
%xlim([0 10]); ylim([0 1])

%% Individual Fits

for i = 1:400;
allBICs(i,1:4) = [d.S1.BSampIR.LN.fits{i}.BICtaeheefix,d.S1.BSampIR.MLN.fits{i}.BICtaeheefix, d.S1.BSampIR.Gam.fits{i}.BICtaeheefix, d.S1.BSampIR.invGam.fits{i}.BICtaeheefix];
end

[~,idx] = min(allBICs, [], 2);

minBICs = zeros(size(allBICs));
minBICs(sub2ind(size(allBICs), (1:size(allBICs,1))', idx)) = 1;

preferDistCounts = sum(minBICs, 1);

aveBICs = mean(allBICs, 1);

%%