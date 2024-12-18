%% REPRODUCING BIGMACS MLN - from BIGMACS histogram data
%Add important paths
addpath('../Functions')

%Load BIGMACS files
lognorm_BIGMACS = readtable("../BIGMACSdata/lognormal.txt");                  % (use x values currently used in BIGMACS)
x = lognorm_BIGMACS.Var1';
BIGMACShist = readmatrix("../BIGMACSdata/Lin2014_sedrateratio_cm_wo_NaN.txt");
BIGMACShist = BIGMACShist(:,4);
TM_BIGMACS = readmatrix("../BIGMACSdata/transition_parameter.txt");

%Note, BIGMACS data is nSR counts, where replicates have already been
%formed to provide weighting by depth

%Set up fitting settings structure
fitS.Lin2014AgeFiltering = 1;
fitS.weighting = "depth"; 
fitS.chi2binN = 10;
fitS.dispChi2 = true;
fitS.mln1RunReps = 1;
fitS.mlnReps = 5;

%Fit Mix Log Norm to data
[MLN_BIGMACS, ~, gmfitBM] = fitMixLogNorm(BIGMACShist, x, 2, 0, fitS.mlnReps);

%Run chi2gof on data vs fitted distribution
[hBMvBM, pBMvBM, chiStatBMvBM] = chi2gof_vsMLN(gmfitBM, log(BIGMACShist), 358, fitS);
if fitS.dispChi2
gcf;
title("BIGMACS data vs BIGMACS mix log normal")
end

figure;
hold on
yyaxis("left")
histogram(BIGMACShist, 'FaceColor', 'k', 'FaceAlpha', 0.1, 'DisplayName', 'BIGMACS nSR counts')
ylabel("Counts")
yyaxis("right")
plot(x, lognorm_BIGMACS.Var2, '-k', 'DisplayName', "BIGMACS", 'LineWidth', 1)
plot(x, MLN_BIGMACS(:,2), '--r', 'DisplayName', "BIGMACS - Reproduction", 'LineWidth', 1)
xlim([0 6])
xlabel("nSR")
ylabel("PDF")
legend()
title("BIGMACS reproduction from weighted nSR counts")
set(gcf, 'Position', [20, 400, 500, 200])

%Can't calculate transition matrix from this data because it does not let
%me know which data comes from which core. Means I don't know when a
%transition actually occurred or when a change in nSR is due to a change in
%core.

%Can't plot the nSR histories of the cores