
%% Fit Mix Log Norm to BIGMACS data for comparison reasons
%Load BIGMACS files
BM.lognorm = readtable("../BIGMACSdata/lognormal.txt");                  % (use x values currently used in BIGMACS)
BM.x = BM.lognorm.Var1';
BM.hist = readmatrix("../BIGMACSdata/Lin2014_sedrateratio_cm_wo_NaN.txt");
BM.hist = BM.hist(:,4);
BM.TM = readmatrix("../BIGMACSdata/transition_parameter.txt");
[MLN_BIGMACS, ~, gmfitBM] = fitMixLogNorm(BM.hist, BM.x, 2, 3);

[BM.fxMLN, BM.pfxMLN] = px_to_pfx(BM.x, BM.lognorm.Var2, @log);

%% Figure of All Data, Newall500 Sampling in nSR, log(nSR) and inv(nSR)
% Load data
dAD = load("../Results/dataT_PFandLin_R200M20_Feb4_fitFeb17.mat");
d = dAD.d;
dAllSampB = dAD.d.Bchron1RunFits.weightedC;
dAllSampN = dAD.d.Newall1RunFits.weightedC;
d2_label = "Newall500";

%Create binedges in invnSR
if max(cell2mat(dAllSampB)) < 15
    binEdges = 0:0.1:15; else
    binEdges = [0:0.1:15, ceil(max(cell2mat(dAllSampB)))];
end

%Create binedges in lognSR
logBinEdges = -5:0.1:5;

%Create binedges in invnSR
if 1/(min(cell2mat(dAllSampB))) < 15
    invBinEdges = 0:0.1:15; else
    invBinEdges = [0:0.1:15, 1./min(cell2mat(dAllSampB))];
end

%How many runs to plot
numruns = 100;

figure;
subplot(3,1,1)
title([d.label ;d2_label])
xlim([0 exp(2)])
hold on
ylabel("cm")
histogram(cell2mat(dAllSampN), 'BinEdges', binEdges)
xlabel("nSR")
subplot(3,1,2)
hold on
histogram(log(cell2mat(dAllSampN)), 'BinEdges', logBinEdges)
xlim([-4 4])
xlabel("log(nSR)")
ylabel("cm")
subplot(3,1,3)
hold on
histogram(1./(cell2mat(dAllSampN)), 'BinEdges', invBinEdges)
xlim([0 (exp(2))])
xlabel("inv(nSR)")
ylabel("cm")
ax = gca;
ax.XDir = 'reverse';


%% Transparent BchronMode, Bchron Samples and Newall500 histograms in log(nSR) space

figure;
subplot(3,1,1)
histogram(log(d.BModeHist), 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts');
xlim([-4 4])
title("BchronMode")
subplot(3,1,2)
for i = 1:numruns
    hold on
    histogram(log(dAllSampB{i}), 'BinEdges', logBinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none')
end
xlim([-4 4])

title("Bchron Samples, Transparent Histograms")
%Now plot each individual histogram with a transparent face color
subplot(3,1,3)
for i = 1:numruns
    hold on
    histogram(log(dAllSampN{i}), 'BinEdges', logBinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none')
end
xlim([-4 4])
xlabel("log(nSR)")
title("Newall 500 Samples, Transparent Histograms")

%% Plot combined histograms in log nSR space
logXLim = [-3 3];
figure;
subplot(3,1,1)
hBchModeLog = histogram(log(d.BModeHist), 'FaceColor', 'k', 'FaceAlpha', 0.2, 'DisplayName', 'Replicated nSR counts');
xlim(logXLim)
ylabel("cm")
title("BchronMode")

subplot(3,1,2)
histogram(log(cell2mat(dAllSampB)), 'BinEdges', logBinEdges)%, 'FaceColor', 'k', 'FaceAlpha', 0.4)
xlim(logXLim)
ylabel("cm")
title("Bchron Samples Combined")

subplot(3,1,3)
histogram(log(cell2mat(dAllSampN)), 'BinEdges', logBinEdges)%, 'FaceColor', 'k', 'FaceAlpha', 0.4)
xlim(logXLim)
xlabel("log(nSR)")
ylabel("cm")
title("Newall 500 Samples Combined")

%% Comparison of new dataset to Lin2014 dataset

dLin = load("../Results/dataT_LinOnly_LinMethod_Dec10_fitFeb17.mat");

%This could be with histograms

%This could simply be with number of 14C dates, cores, and cm
disp("Lin data has " + num2str(length(dLin.d.dataT.cores)) + " cores, with " + num2str(sum(dLin.d.dataT.num14cpairs)) + " age pairs, with a total sediment length of " + num2str(sum(dLin.d.dataT.sedimentlength)) + "cm")
disp("Chosen data has " + num2str(length(dAD.d.dataT.cores)) + " cores, with " + num2str(sum(dAD.d.dataT.num14cpairs)) + " age pairs, with a total sediment length of " + num2str(sum(dAD.d.dataT.sedimentlength)) + "cm")
%Need to mention that not all cores from Lin et al., 2014, are being used
%in my new dataset

%% Chi-squared GOF test on 1 sampling of Newall500 sampling approach - with evenly spaced bars (except for end bars)
rng(1)
r = randi(numruns, 1);

f_log = @(x) log(x);
OneNewSamp = dAD.d.Newall1RunFits.weightedC{r};
OneNewMLN = dAD.d.Newall1RunFits.MLNPDFs(:,1);
OneNewMLNpdf.numParams = 6;
[OneNewMLNpdf.x, OneNewMLNpdf.px] = px_to_pfx(dAD.d.x,OneNewMLN, f_log);
OneNewInvGam = dAD.d.Newall1RunFits.invGamPDFs(:,1);
OneNewInvGampdf.numParams = 2;
[OneNewInvGampdf.x, OneNewInvGampdf.px] = px_to_pfx(dAD.d.x,OneNewInvGam, f_log);
fitS.dispChi2 = true;

[h1, p1, chiStat1, h2, p2, chiStat2] = chi2_dataVStwopdfVECs(log(dAD.d.Newall1RunFits.weightedC{r}), dAD.d.Newall1RunFits.numCpairs(r), 20,OneNewMLNpdf , OneNewInvGampdf, fitS);

% Need to make histogram bars equally spaced

%% Compare newall sampled data and fits to Bchron sampled data and fits

%% Check fits with changing restriction on Newall sampling approach

%% Compare histograms from depth-weighted, age-weighted, unweighted

%% Check impact of using single mean SR vs actual mean SR for sampling

