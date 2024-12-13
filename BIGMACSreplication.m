%% Assess my attempts of replicating BIGMACS
%Add important paths
addpath('Functions')

%Load BIGMACS files
lognorm_BIGMACS = readtable("BIGMACSdata/lognormal.txt");                  % (use x values currently used in BIGMACS)
x = lognorm_BIGMACS.Var1';
BIGMACShist = readmatrix("BIGMACSdata/Lin2014_sedrateratio_cm_wo_NaN.txt");
BIGMACShist = BIGMACShist(:,4);
TM_BIGMACS = readmatrix("BIGMACSdata/transition_parameter.txt");

%Note, BIGMACS data is nSR counts, where replicates have already been
%formed to provide weighting by depth

%Set up fitting settings structure
fitS.Lin2014AgeFiltering = 1;
fitS.weighting = "depth"; 
fitS.chi2binN = 10;
fitS.dispChi2 = false;

%Fit Mix Log Norm to data
[MLN_BIGMACS, ~, gmfitBM] = fitMixLogNorm(BIGMACShist, x, 2, 0, 5);

%Run chi2gof on data vs fitted distribution
disp("BIGMACS data vs BIGMACS mix log normal");
[hBMvBM, pBMvBM, chiStatBMvBM] = chi2gof_vsMLN(gmfitBM, log(BIGMACShist), 358, fitS);
if fitS.dispChi2
gcf;
title("BIGMACS data vs BIGMACS mix log normal")
end

%Can't calculate transition matrix from this data because it does not let
%me know which data comes from which core. Means I don't know when a
%transition actually occurred or when a change in nSR is due to a change in
%core.

%Can't plot the nSR histories of the cores

%% Use my results using Lin 2014 cores and replicated methodology
%Load my results file
load("Results/dataT_LinOnly_LinMethod_Dec10.mat")
LinCoreLinMeth.dataT = dataT;
LinCoreLinMeth.S = S;

sizeCores = numel(LinCoreLinMeth.dataT.lats);
numCores = sizeCores;

%Fit MLN distribution to my Bchron Mode data
disp("LinMethod (mode) data vs best fit mln")
[LinCoreLinMeth.mixLogBMode, LinCoreLinMeth.BModeHist,~,~,~,LinCoreLinMeth.gmfitBmode, LinCoreLinMeth.ncBmode, LinCoreLinMeth.h, LinCoreLinMeth.p,LinCoreLinMeth.chiStat] = plotSRandResHistograms(dataT.bchronMode, x, true(sizeCores), 3, 1, 2, 0, "", 1, fitS);
if fitS.dispChi2
gcf; title("chi2gof: LinCoreBmode Data vs Best Fit MLN")
end
%See how my distributions perform with chi2gof on BIGMACS gmfit
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

%% Try getting chi2gof for the individual Bchron runs
numruns = 400;
fitS.dispChi2 = false;
[LinCoreLinMeth.MLN1R.pdfs, LinCoreLinMeth.MLN1R.c95up, LinCoreLinMeth.MLN1R.c95down, LinCoreLinMeth.MLN1R.mus, LinCoreLinMeth.MLN1R.sigmas, LinCoreLinMeth.MLN1R.outputS] = SingleRunLogNorms(LinCoreLinMeth.dataT.bchronProb, true(sizeCores, 1), numruns, x, 2, 3, 4, 0, fitS);

% Test the chi2gof of each fit to BIGMACS
    h1R = NaN(numruns,1);
    p1R = NaN(numruns,1);
    fitS.dispChi2 = false;
    chiStat1RunT = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);

for i = 1:numruns
    [h1R(i), p1R(i), chiStat1Run] = chi2gof_vsMLN(gmfitBM, log(LinCoreLinMeth.MLN1R.outputS.weightedC{i}), LinCoreLinMeth.MLN1R.outputS.numCpairs(i), fitS);
    chiStat1RunT = [chiStat1RunT; struct2table(chiStat1Run)]; %#ok<AGROW>
end

chiStat1RunT = addvars(chiStat1RunT, h1R, p1R, 'Before', "chi2stat");

LinCoreLinMeth.MLN1R.chiStat1RunT = chiStat1RunT;

%% Bring in results using Bchron but with updated core dataset

load("Results/dataT_LinandPF_LinMethod_Dec10.mat")

sizeCores = numel(dataT.lats);
numCores = sizeCores;

%Fit MLN distribution to my Bchron Mode data
disp("LinMethod (mode) data vs best fit mln")
[NewCoreLinMeth.mixLogBMode, NewCoreLinMeth.BModeHist,~,~,~,NewCoreLinMeth.gmfitBmode, NewCoreLinMeth.ncBmode, NewCoreLinMeth.h, NewCoreLinMeth.p, NewCoreLinMeth.chiStat] = plotSRandResHistograms(dataT.bchronMode, x, true(sizeCores), 3, 1, 2, 0, "", 1, fitS);
if fitS.dispChi2
gcf; title("chi2gof: NewCoreBmode Data vs Best Fit MLN")
end
%See how my distributions perform with chi2gof on BIGMACS gmfit
disp("LinMethod data (mode) vs BIGMACS MLN")
[NewCoreLinMeth.hBM, NewCoreLinMeth.pBM, NewCoreLinMeth.chiStatBM] = chi2gof_vsMLN(gmfitBM, log(NewCoreLinMeth.BModeHist), NewCoreLinMeth.ncBmode, fitS);
if fitS.dispChi2
gcf; title("chi2gof: NewCoreBmode Data vs BIGMACS fit")
end
%Calculate TM of my Bchron Mode data
S.weighting = "age";
[~,~,NewCoreLinMeth.TMnewBchNoWeight, NewCoreLinMeth.TMnewBchWeightD] = TMcalculation(dataT.bchronMode, true(sizeCores), S);

orderInd = 1:numCores;
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

%% Try getting chi2gof for the individual Bchron runs
numruns = 400;
fitS.dispChi2 = false;
[NewCoreLinMeth.MLN1R.pdfs, NewCoreLinMeth.MLN1R.c95up, NewCoreLinMeth.MLN1R.c95down, NewCoreLinMeth.MLN1R.mus, NewCoreLinMeth.MLN1R.sigmas, NewCoreLinMeth.MLN1R.outputS] = SingleRunLogNorms(NewCoreLinMeth.dataT.bchronProb, true(sizeCores, 1), numruns, x, 2, 3, 4, 0, fitS);

% Test the chi2gof of each fit to BIGMACS
    h1R = NaN(1,numruns);
    p1R = NaN(1,numruns);
    fitS.dispChi2 = false;
    chiStat1RunT = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);

for i = 1:numruns
    [h1R(i), p1R(i), chiStat1Run] = chi2gof_vsMLN(gmfitBM, log(NewCoreLinMeth.MLN1R.outputS.weightedC{i}), NewCoreLinMeth.MLN1R.outputS.numCpairs(i), fitS);
    chiStat1RunT = [chiStat1RunT; struct2table(chiStat1Run)]; %#ok<AGROW>
end

chiStat1RunT = addVars(chiStat1RunT, h1R, p1R, 'Before', "chi2stat");

NewCoreLinMeth.MLN1R.chiStat1RunT = chiStat1RunT;

%% Create table to show chi2Stats for each pairing
numPDFs = 3;
numDatas = 3;
chi2Table = table('Size', [numPDFs, numDatas], 'VariableTypes', ["double", "double", "double"],'VariableNames',["BMhist", "LinCoresLinMethBMode", "NewCoresLinMethBMode"], 'RowNames',["BMpdf", "LCLMBModeMLN", "NCLMBModeMLN"]);
chi2Table.BMhist(1) = chiStatBMvBM.chi2stat;
chi2Table.LinCoresLinMethBMode(1) = LinCoreLinMeth.chiStatBM.chi2stat;
chi2Table.LinCoresLinMethBMode(2) = LinCoreLinMeth.chiStat.chi2stat;
chi2Table.NewCoresLinMethBMode(1) = NewCoreLinMeth.chiStatBM.chi2stat;
chi2Table.NewCoresLinMethBMode(3) = NewCoreLinMeth.chiStat.chi2stat;


