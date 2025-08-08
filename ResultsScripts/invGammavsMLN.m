%% Bring in data with new dataset but with Marine 20 and 200 year R uncertainty
%Add important paths
addpath('../Functions')

%Load BIGMACS files
BM.lognorm = readtable("../BIGMACSdata/lognormal.txt");                  % (use x values currently used in BIGMACS)
BM.x = BM.lognorm.Var1';
BM.hist = readmatrix("../BIGMACSdata/Lin2014_sedrateratio_cm_wo_NaN.txt");
BM.hist = BM.hist(:,4);
BM.TM = readmatrix("../BIGMACSdata/transition_parameter.txt");

%Load results file
load("../Results/dataT_PFandLin_R200M20_Feb4.mat")
d.dataT = dataT;
d.S = S;
d.numCores = numel(d.dataT.lats);

%Set up fitS structure
fitS.dispChi2 = false;
fitS.Lin2014AgeFiltering = true;
fitS.mlnReps = 5;
fitS.weighting = "depth";
fitS.chi2binN = 200;
fitS.invXbinEdges = 0:0.1:15;
fitS.enforceBinSizeLimits = true;

% Want to fit some inverse gamma distributions and mix log normal
% distributions to the same data to see how they look, and compare their
% chi2stat values...

%Fit MLN distribution to my Bchron Mode data
disp("data vs best fit mln")
[d.mixLogBMode, d.BModeHist,~,~,~,...
    d.gmfitBmode, d.ncBmode, d.h,...
    d.p, d.chiStat]...
    = plotSRandResHistograms(d.dataT.bchronMode, BM.x, true(d.numCores), 3, 1,...
    2, 0, "", 0, fitS);
if fitS.dispChi2
    gcf; title("chi2gof: Data vs Best Fit MLN")
end

%Fit Inverse gamma distribution to my Bchron Mode data
[d.invGamBMode, d.invGamBModeProb] = fitGamma2invSR(d.dataT.bchronMode, true(d.numCores), fitS);

%Plot both inverse gamma and mln
figure()
hold on
yyaxis left
plot(d.invGamBMode, d.invGamBModeProb, '--r', 'LineWidth', 1, 'DisplayName', 'MLN');
plot(d.mixLogBMode(:,1), d.mixLogBMode(:,2), '--k', 'LineWidth', 1, 'DisplayName', 'Inverse Gamma');
ylabel("pdf")
yyaxis right
histogram(d.BModeHist, 'FaceColor','b', 'FaceAlpha', 0.1, 'DisplayName', 'Sediment (cm)')
ylabel("cm")
xlabel('nSR')
xlim([0 6])
legend()

%Plot both on log(nSR) space
%Need to change inverse gamma pdf to log(nSR) space
f_log = @(x) log(x);
[d.invGamBMode_lnSR, d.invGamBModeProb_lnSR] = px_to_pfx(d.invGamBMode, d.invGamBModeProb, f_log);
[d.MLNBMode_lnSR, d.MLNBModeProb_lnSR] = px_to_pfx(d.mixLogBMode(:,1), d.mixLogBMode(:,2), f_log);

figure()
hold on
yyaxis left
plot(d.MLNBMode_lnSR, d.MLNBModeProb_lnSR, '-k', 'LineWidth', 1, 'DisplayName', 'MLN');
plot(d.invGamBMode_lnSR, d.invGamBModeProb_lnSR, '-r', 'LineWidth', 1, 'DisplayName', 'Inverse Gamma')
ylabel("pdf")
yyaxis right
histogram(log(d.BModeHist), 'FaceColor','b', 'FaceAlpha', 0.1, 'DisplayName', 'Sediment (cm)')
ylabel("cm")
xlabel('log(nSR)')
legend()

%Plot both on log(nSR) space
%Need to change inverse gamma pdf to log(nSR) space
f_inv = @(x) 1./x;
[d.invGamBMode_invnSR, d.invGamBModeProb_invnSR] = px_to_pfx(d.invGamBMode, d.invGamBModeProb, f_inv);
[d.MLNBMode_invnSR, d.MLNBModeProb_invnSR] = px_to_pfx(d.mixLogBMode(:,1), d.mixLogBMode(:,2), f_inv);

figure()
hold on
yyaxis left
plot(d.MLNBMode_invnSR, d.MLNBModeProb_invnSR, '-k', 'LineWidth', 1, 'DisplayName', 'MLN');
plot(d.invGamBMode_invnSR, d.invGamBModeProb_invnSR, '-r', 'LineWidth', 1, 'DisplayName', 'Inverse Gamma')
ylabel("pdf")
yyaxis right
histogram(1./d.BModeHist, 'FaceColor','b', 'FaceAlpha', 0.1, 'DisplayName', 'Sediment (cm)')
ylabel("cm")
xlabel('inverse(nSR)')
legend()
xlim([0 6])

%% Chi squared using the pdf vector of MLN and InvGam on log nsR
d.pdfVEC1.x = d.MLNBMode_lnSR;
d.pdfVEC1.fx = d.MLNBModeProb_lnSR;
d.pdfVEC1.numParams = 6;
fitS.dispChi2 = true;
d.pdfVEC2.x = d.invGamBMode_lnSR;
d.pdfVEC2.fx = d.invGamBModeProb_lnSR;
d.pdfVEC2.numParams = 2;
chi2_dataVStwopdfVECs(log(d.BModeHist), d.ncBmode, 20, d.pdfVEC1, d.pdfVEC2, fitS)

%%

fitS.OneRun.weightRepDP = 3;
fitS.OneRun.weightRepInflator = 4;
fitS.OneRun.MLNReps = 2;
fitS.dispChi2 = false;
[outS] = SRun_MLNandInvGam(d.dataT.bchronProb, true(d.numCores,1), 100, BM.x, 2, fitS);

