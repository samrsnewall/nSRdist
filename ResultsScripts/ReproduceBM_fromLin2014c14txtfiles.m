%% REPRODUCING BIGMACS MLN - Starting from Lin et al radiocarbon files
%Add important paths
addpath('../Functions')

%Load results file
load("../Results/dataT_LinOnly_LinMethod_Dec10.mat")
%Store variables in structure
LinCoreLinMeth.dataT = dataT;
LinCoreLinMeth.S = S;

%Load BIGMACS files
lognorm_BIGMACS = readtable("../BIGMACSdata/lognormal.txt");                  % (use x values currently used in BIGMACS)
x = lognorm_BIGMACS.Var1';
BIGMACShist = readmatrix("../BIGMACSdata/Lin2014_sedrateratio_cm_wo_NaN.txt");
BIGMACShist = BIGMACShist(:,4);
TM_BIGMACS = readmatrix("../BIGMACSdata/transition_parameter.txt");



%Find number of cores
sizeCores = numel(LinCoreLinMeth.dataT.lats);
numCores = sizeCores;

%Fit MLN distribution to my Bchron Mode data
disp("LinMethod (mode) data vs best fit mln")
[LinCoreLinMeth.mixLogBMode, LinCoreLinMeth.BModeHist,~,~,~,LinCoreLinMeth.gmfitBmode, LinCoreLinMeth.ncBmode, LinCoreLinMeth.h, LinCoreLinMeth.p,LinCoreLinMeth.chiStat] = plotSRandResHistograms(dataT.bchronMode, x, true(sizeCores), 3, 1, 2, 0, "", 1, fitS);
if fitS.dispChi2
gcf; title("chi2gof: LinCoreBmode Data vs Best Fit MLN")
end

%Plot new histogram and MLN alongside BIGMACS MLN
figure;
hold on
yyaxis("left")
histogram(LinCoreLinMeth.BModeHist, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Reproduced nSR counts')
ylabel("Counts")
yyaxis("right")
plot(x, lognorm_BIGMACS.Var2, '-k', 'DisplayName', "BIGMACS", 'LineWidth', 1)
plot(x,LinCoreLinMeth.mixLogBMode(:,2), '-r', 'DisplayName', "Reproduced MLN", 'LineWidth', 1)
xlim([0 6])
xlabel("nSR")
ylabel("PDF")
legend()
title("BIGMACS reproduction from C14 files")
set(gcf, 'Position', [20, 400, 500, 200])

%See how my distribution performs with chi2gof on BIGMACS gmfit
disp("LinMethod data (mode) vs BIGMACS MLN")
[LinCoreLinMeth.hBM, LinCoreLinMeth.pBM, LinCoreLinMeth.chiStatBM] = chi2gof_vsMLN(gmfitBM, log(LinCoreLinMeth.BModeHist), LinCoreLinMeth.ncBmode, fitS);
if fitS.dispChi2
gcf; title("chi2gof: LinCoreBmode Data vs BIGMACS fit")
end

%Calculate TM of my Bchron Mode data
S.weighting = "age";
[~,~,LinCoreLinMeth.TMnewBchNoWeight, LinCoreLinMeth.TMnewBchWeightD] = TMcalculation(dataT.bchronMode, true(sizeCores), S);

%% Plot nSR histories
%Plot nSR histories of all the cores (Recreating plot S4 from Lin et al.,
%2014)

lin2014figS4Order = ["MD84-527", "MD88-770", "SO42-74KL", "GeoB7920-2", "GeoB9508-5", "GeoB9526", "KF13", "KNR31-GPC5", "MD03-2698", "MD99-2334", "SU81-18", "H214", "MD01-2421", "ODP1145", "SO50-31KL", "TR163-22", "V19-30", "W8709A-8", "M35003-4", "MD99-2339", "POS200_10_6-2", "DSDP594", "GIK17961-2", "GIK17964-2", "MD97-2120","MD97-2151","GeoB1711-4", "KNR159-5-36GGC", "MD95-2042","MD01-2416", "MD02-2489", "MD98-2181", "W8709A-13", "MD07-3076", "RC11-83", "GIK17940-2", "V35-5"];
[~, orderInd] = ismember(lin2014figS4Order, dataT.cores);
figure;
for i = 1:numCores
    subplot(ceil(numCores./5), 5, i)
    nSRs = dataT.bchronMode{orderInd(i)}(1,2:end);
    ages = cumsum(dataT.bchronMode{orderInd(i)}(4,:))./1000;
    stairs(ages, [nSRs, nSRs(end)], '-b')
    set(gca, 'YScale', 'log')
    ylim([0.1, 10])
    xlim([0 45])
    title(dataT.cores(orderInd(i)))
end
fontsize(gcf, "scale", 0.6)

%From looking at my first attempt of the core nSR history figure, I realized that the mode
%of the Bchron ages are not all motonically increasing - there was a
%negative nSR returned by my BchronMode method in V35-5 and MD02-2489.
%However, I also noticed that the stair plot in figure S4 of Lin et al.,
%2014 isn't actually a strict stair plot. There are lines which have
%non-zero gradients, and some that show some odd behaviour, such as
%MD95-2042 between ages 25-35kyr, MD98-2181 at age 12kyr and MD84-527 at
%10kyr. Perhaps a result of creating the stair plot after filtering for the
%age gaps (which I did not do before plotting the histories). The
%sedimentation rate data they use is actually provided for each of the
%cores (in a useless and frustrating format). However, having looked at
%that, and how the SR values of some (but not all) closely spaced dates are
%the exact same, perhaps they are using a single Bchron run, instead of the
%Bchron mode of the ages... This would make sense because each Bchron run
%is a set of SR rates applied for a varying amount of time, so they have
%areas of consistent SR across multiple dates (which is not the case for
%the mode of the Bchronology results)!

