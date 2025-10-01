
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

%% Figure of All Data, Newall500 Sampling in nSR, log(nSR) and inv(nSR)
% Load data
dAD = load("../Results/dataT_PFandLin_R200M20_Feb4_fitFeb18.mat");
d = dAD.d;
dN500_label = "Newall500";

%Create binedges in invnSR
if max(cell2mat(dAD.d.S1.BChIR.weightedC)) < 15
    binEdges = 0:0.1:15; else
    binEdges = [0:0.1:15, ceil(max(cell2mat(dAD.d.S1.BChIR.weightedC)))];
end

%Create binedges in lognSR
logBinEdges = -5:0.1:5;

%Create binedges in invnSR
if 1/(min(cell2mat(dAD.d.S1.BChIR.weightedC))) < 15
    invBinEdges = 0:0.1:15; else
    invBinEdges = [0:0.1:15, 1./min(cell2mat(dAD.d.S1.BChIR.weightedC))];
end

%How many runs to plot
numruns = 100;

figure;
subplot(3,1,1)
title([d.label ;dN500_label])
xlim([0 exp(2)])
hold on
ylabel("cm")
histogram(cell2mat(dAD.d.S1.New500IR.weightedC), 'BinEdges', binEdges)
xlabel("nSR")
subplot(3,1,2)
hold on
histogram(log(cell2mat(dAD.d.S1.New500IR.weightedC)), 'BinEdges', logBinEdges)
xlim([-4 4])
xlabel("log(nSR)")
ylabel("cm")
subplot(3,1,3)
hold on
histogram(1./(cell2mat(dAD.d.S1.New500IR.weightedC)), 'BinEdges', invBinEdges)
xlim([0 (exp(2))])
xlabel("inv(nSR)")
ylabel("cm")
ax = gca;
ax.XDir = 'reverse';


%% Transparent BchronMode, Bchron Samples and Newall500 histograms in log(nSR) space

figure;
subplot(3,1,1)
histogram(log(dAD.d.S1.BMode.Hist), 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts');
xlim([-4 4])
title("BchronMode")
subplot(3,1,2)
for i = 1:numruns
    hold on
    histogram(log(dAD.d.S1.BChIR.weightedC{i}), 'BinEdges', logBinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none')
end
xlim([-4 4])

title("Bchron Samples, Transparent Histograms")
%Now plot each individual histogram with a transparent face color
subplot(3,1,3)
for i = 1:numruns
    hold on
    histogram(log(dAD.d.S1.New500IR.weightedC{i}), 'BinEdges', logBinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none')
end
xlim([-4 4])
xlabel("log(nSR)")
title("Newall 500 Samples, Transparent Histograms")

%% Plot combined histograms in log nSR space
logXLim = [-3 3];
figure;
subplot(3,1,1)
hBchModeLog = histogram(log(dAD.d.S1.BMode.Hist), 'FaceColor', 'k', 'FaceAlpha', 0.2, 'DisplayName', 'Replicated nSR counts');
xlim(logXLim)
ylabel("cm")
title("BchronMode")

subplot(3,1,2)
histogram(log(cell2mat(dAD.d.S1.BChIR.weightedC)), 'BinEdges', logBinEdges)%, 'FaceColor', 'k', 'FaceAlpha', 0.4)
xlim(logXLim)
ylabel("cm")
title("Bchron Samples Combined")

subplot(3,1,3)
histogram(log(cell2mat(dAD.d.S1.New500IR.weightedC)), 'BinEdges', logBinEdges)%, 'FaceColor', 'k', 'FaceAlpha', 0.4)
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
rng(4)
r = randi(numruns, 1);

f_log = @(x) log(x);
OneNewSamp = dAD.d.S1.New500IR.weightedC{r};
OneNewMLN = dAD.d.S1.New500IR.MLNPDFs(:,1);
OneNewMLNpdf.numParams = 6;
[OneNewMLNpdf.x, OneNewMLNpdf.px] = px_to_pfx(dAD.d.x,OneNewMLN, f_log);
OneNewInvGam = dAD.d.S1.New500IR.invGamPDFs(:,1);
OneNewInvGampdf.numParams = 2;
[OneNewInvGampdf.x, OneNewInvGampdf.px] = px_to_pfx(dAD.d.x,OneNewInvGam, f_log);
fitS.dispChi2 = true;

[h1, p1, chiStat1, h2, p2, chiStat2] = chi2_dataVStwopdfVECs(log(dAD.d.S1.New500IR.weightedC{r}), dAD.d.S1.New500IR.numCpairs(r), 20,OneNewMLNpdf , OneNewInvGampdf, fitS);

%% Compare newall sampled data and fits to Bchron sampled data and fits

datasetB1 = dAD.d.S1.BChIR.weightedC;
dB1_label = "Bchron Sampling";

datasetN500 = dAD.d.S1.New500IR.weightedC;
dN500_label = "Newall500";

dAD.d.S1.BChIR.MLN.lnSR.x = NaN(length(d.x), numruns);
dAD.d.S1.BChIR.MLN.lnSR.px = NaN(length(d.x), numruns);
for i = 1:100
    [dAD.d.S1.BChIR.MLN.lnSR.x(:,i), dAD.d.S1.BChIR.MLN.lnSR.px(:,i)] = px_to_pfx(d.x, dAD.d.S1.BChIR.MLNPDFs(:,i), @log);
end

dAD.d.S1.BChIR.invGam.lnSR.x = NaN(length(d.x), numruns);
dAD.d.S1.BChIR.invGam.lnSR.px = NaN(length(d.x), numruns);
for i = 1:100
    [dAD.d.S1.BChIR.invGam.lnSR.x(:,i), dAD.d.S1.BChIR.invGam.lnSR.px(:,i)] = px_to_pfx(d.x, dAD.d.S1.BChIR.invGamPDFs(:,i), @log);
end

%Find which distribution fits better for each run
chooseMLN_B1 = sum(dAD.d.S1.BChIR.MLNchiStats.p > dAD.d.S1.BChIR.InvGamChiStats.p);
chooseInvGam_B1 = numruns - chooseMLN_B1;

dAD.d.S1.New500IR.MLN.lnSR.x = NaN(length(d.x), numruns);
dAD.d.S1.New500IR.MLN.lnSR.px = NaN(length(d.x), numruns);
for i = 1:100
    [dAD.d.S1.New500IR.MLN.lnSR.x(:,i), dAD.d.S1.New500IR.MLN.lnSR.px(:,i)] = px_to_pfx(d.x, dAD.d.S1.BChIR.MLNPDFs(:,i), @log);
end

dAD.d.S1.New500IR.invGam.lnSR.x = NaN(length(d.x), numruns);
dAD.d.S1.New500IR.invGam.lnSR.px = NaN(length(d.x), numruns);
for i = 1:100
    [dAD.d.S1.New500IR.invGam.lnSR.x(:,i), dAD.d.S1.New500IR.invGam.lnSR.px(:,i)] = px_to_pfx(d.x, dAD.d.S1.BChIR.invGamPDFs(:,i), @log);
end

%Find which distribution fits better for each run
chooseMLN_Newall500 = sum(dAD.d.S1.New500IR.MLNchiStats.p > dAD.d.S1.New500IR.InvGamChiStats.p);
chooseInvGam_Newall500 = numruns - chooseMLN_Newall500;

figure()
subplot(2,1,1)
xlim([-3 3])
yyaxis left
hold on
ylabel("cm")
for i = 1:numruns
    hold on
    if i ==1 
    histogram(log(dAD.d.S1.BChIR.weightedC{i}), 'BinEdges', logBinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'on', 'DisplayName', 'Histograms')
    else
    histogram(log(dAD.d.S1.BChIR.weightedC{i}), 'BinEdges', logBinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'off')
    end   
end
yyaxis right
hold on
plot(dAD.d.S1.BChIR.MLN.lnSR.x, dAD.d.S1.BChIR.MLN.lnSR.px,'LineStyle', '-', 'Color', [0,0,1,0.1],'Marker','none' ,'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [0,0,1], 'Marker','none', 'DisplayName', "MLNs: " + num2str(chooseMLN_B1))
plot(dAD.d.S1.BChIR.invGam.lnSR.x, dAD.d.S1.BChIR.invGam.lnSR.px,'LineStyle', '-', 'Color', [1,0,0,0.1],'Marker','none', 'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [1,0,0], 'Marker','none', 'DisplayName', "InvGammas: " + num2str(chooseInvGam_B1))
ylim([0 1])
xlabel("nSR")
ylabel("PDF")
legend()
title("Bchron Samplings")
subplot(2,1,2)
xlim([-3 3])
yyaxis left
hold on
ylabel("cm")
for i = 1:numruns
    hold on
    if i ==1 
    histogram(log(dAD.d.S1.New500IR.weightedC{i}), 'BinEdges', logBinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'on', 'DisplayName', 'Histograms')
    else
    histogram(log(dAD.d.S1.New500IR.weightedC{i}), 'BinEdges', logBinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'off')
    end   
end
yyaxis right
hold on
plot(dAD.d.S1.New500IR.MLN.lnSR.x, dAD.d.S1.New500IR.MLN.lnSR.px,'LineStyle', '-', 'Color', [0,0,1,0.1],'Marker','none' ,'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [0,0,1], 'Marker','none', 'DisplayName', "MLNs: " + num2str(chooseMLN_Newall500))
plot(dAD.d.S1.New500IR.invGam.lnSR.x, dAD.d.S1.New500IR.invGam.lnSR.px,'LineStyle', '-', 'Color', [1,0,0,0.1],'Marker','none', 'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [1,0,0], 'Marker','none', 'DisplayName', "InvGammas: " + num2str(chooseInvGam_Newall500))
ylim([0 1])
xlabel("nSR")
ylabel("PDF")
title("Newall Samplings")
legend()

%% Check fits with changing restriction on Newall sampling approach

%Convert all Newall0 pdfs to log space
datasetN0 = dAD.d.S1.New0IR.weightedC;
dN0_label = "Newall0";

dAD.d.S1.New0IR.MLN.lnSR.x = NaN(length(d.x), numruns);
dAD.d.S1.New0IR.MLN.lnSR.px = NaN(length(d.x), numruns);
for i = 1:numruns
    [dAD.d.S1.New0IR.MLN.lnSR.x(:,i), dAD.d.S1.New0IR.MLN.lnSR.px(:,i)] = px_to_pfx(d.x, dAD.d.S1.New0IR.MLNPDFs(:,i), @log);
end

dAD.d.S1.New0IR.invGam.lnSR.x = NaN(length(d.x), numruns);
dAD.d.S1.New0IR.invGam.lnSR.px = NaN(length(d.x), numruns);
for i = 1:numruns
    [dAD.d.S1.New0IR.invGam.lnSR.x(:,i), dAD.d.S1.New0IR.invGam.lnSR.px(:,i)] = px_to_pfx(d.x, dAD.d.S1.New0IR.invGamPDFs(:,i), @log);
end

%Find which distribution fits better for each run
chooseMLN_Newall0 = sum(dAD.d.S1.New0IR.MLNchiStats.p > dAD.d.S1.New0IR.InvGamChiStats.p);
chooseInvGam_Newall0 = numruns - chooseMLN_Newall0;

%Convert all Newall1000 pdfs to log space
dAD.d.S1.New1000IR.MLN.lnSR.x = NaN(length(d.x), numruns);
dAD.d.S1.New1000IR.MLN.lnSR.px = NaN(length(d.x), numruns);
for i = 1:numruns
    [dAD.d.S1.New1000IR.MLN.lnSR.x(:,i), dAD.d.S1.New1000IR.MLN.lnSR.px(:,i)] = px_to_pfx(d.x, dAD.d.S1.New1000IR.MLNPDFs(:,i), @log);
end

dAD.d.S1.New1000IR.invGam.lnSR.x = NaN(length(d.x), numruns);
dAD.d.S1.New1000IR.invGam.lnSR.px = NaN(length(d.x), numruns);
for i = 1:numruns
    [dAD.d.S1.New1000IR.invGam.lnSR.x(:,i), dAD.d.S1.New1000IR.invGam.lnSR.px(:,i)] = px_to_pfx(d.x, dAD.d.S1.New1000IR.invGamPDFs(:,i), @log);
end

%Find which distribution fits better for each run
chooseMLN_Newall1000 = sum(dAD.d.S1.New1000IR.MLNchiStats.p > dAD.d.S1.New1000IR.InvGamChiStats.p);
chooseInvGam_Newall1000 = numruns - chooseMLN_Newall1000;


figure()
subplot(3,1,1)
xlim([-3 3])
yyaxis left
hold on
ylabel("cm")
for i = 1:numruns
    hold on
    if i ==1 
    histogram(log(dAD.d.S1.New0IR.weightedC{i}), 'BinEdges', logBinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'on', 'DisplayName', 'Histograms')
    else
    histogram(log(dAD.d.S1.New0IR.weightedC{i}), 'BinEdges', logBinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'off')
    end   
end
yyaxis right
hold on
plot(dAD.d.S1.New0IR.MLN.lnSR.x, dAD.d.S1.New0IR.MLN.lnSR.px,'LineStyle', '-', 'Color', [0,0,1,0.1],'Marker','none' ,'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [0,0,1], 'Marker','none', 'DisplayName', "MLNs: " + num2str(chooseMLN_Newall0))
plot(dAD.d.S1.New0IR.invGam.lnSR.x, dAD.d.S1.New0IR.invGam.lnSR.px,'LineStyle', '-', 'Color', [1,0,0,0.1],'Marker','none', 'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [1,0,0], 'Marker','none', 'DisplayName', "InvGammas: " + num2str(chooseInvGam_Newall0))
ylim([0 1])
xlabel("nSR")
ylabel("PDF")
legend()
title("Newall0 Samplings")
subplot(3,1,2)
xlim([-3 3])
yyaxis left
hold on
ylabel("cm")
for i = 1:numruns
    hold on
    if i ==1 
    histogram(log(dAD.d.S1.New500IR.weightedC{i}), 'BinEdges', logBinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'on', 'DisplayName', 'Histograms')
    else
    histogram(log(dAD.d.S1.New500IR.weightedC{i}), 'BinEdges', logBinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'off')
    end   
end
yyaxis right
hold on
plot(dAD.d.S1.New500IR.MLN.lnSR.x, dAD.d.S1.New500IR.MLN.lnSR.px,'LineStyle', '-', 'Color', [0,0,1,0.1],'Marker','none' ,'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [0,0,1], 'Marker','none', 'DisplayName', "MLNs: " + num2str(chooseMLN_Newall500))
plot(dAD.d.S1.New500IR.invGam.lnSR.x, dAD.d.S1.New500IR.invGam.lnSR.px,'LineStyle', '-', 'Color', [1,0,0,0.1],'Marker','none', 'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [1,0,0], 'Marker','none', 'DisplayName', "InvGammas: " + num2str(chooseInvGam_Newall500))
ylim([0 1])
xlabel("nSR")
ylabel("PDF")
title("Newall500 Samplings")
legend()
subplot(3,1,3)
xlim([-3 3])
yyaxis left
hold on
ylabel("cm")
for i = 1:numruns
    hold on
    if i ==1 
    histogram(log(dAD.d.S1.New1000IR.weightedC{i}), 'BinEdges', logBinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'on', 'DisplayName', 'Histograms')
    else
    histogram(log(dAD.d.S1.New1000IR.weightedC{i}), 'BinEdges', logBinEdges, 'FaceAlpha', 0.05, 'FaceColor', 'k', 'EdgeColor','none', 'HandleVisibility', 'off')
    end   
end
yyaxis right
hold on
plot(dAD.d.S1.New1000IR.MLN.lnSR.x, dAD.d.S1.New1000IR.MLN.lnSR.px,'LineStyle', '-', 'Color', [0,0,1,0.1],'Marker','none' ,'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [0,0,1], 'Marker','none', 'DisplayName', "MLNs: " + num2str(chooseMLN_Newall500))
plot(dAD.d.S1.New1000IR.invGam.lnSR.x, dAD.d.S1.New1000IR.invGam.lnSR.px,'LineStyle', '-', 'Color', [1,0,0,0.1],'Marker','none', 'HandleVisibility', 'off')
plot(NaN, NaN,'LineStyle', '-', 'Color', [1,0,0], 'Marker','none', 'DisplayName', "InvGammas: " + num2str(chooseInvGam_Newall500))
ylim([0 1])
xlabel("nSR")
ylabel("PDF")
title("Newall1000 Samplings")
legend()

%% Check impact of using single mean SR vs actual mean SR for sampling


%% Compare histograms from depth-weighted, age-weighted, unweighted


