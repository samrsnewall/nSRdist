
%Add important paths
addpath('Functions')

%% REPRODUCING BIGMACS MLN - from BIGMACS histogram data

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
fitS.dispChi2 = true;
fitS.mln1RunReps = 1;
fitS.mlnReps = 5;

%Fit Mix Log Norm to data
[MLN_BIGMACS, ~, gmfitBM] = fitMixLogNorm(BIGMACShist, x, 2, 0, 5);

%Run chi2gof on data vs fitted distribution
disp("BIGMACS data vs BIGMACS mix log normal");
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


%% REPRODUCING BIGMACS MLN - Starting from Lin et al radiocarbon files
%Load correct results file
load("Results/dataT_LinOnly_LinMethod_Dec10.mat")

%Store variables in structure
LinCoreLinMeth.dataT = dataT;
LinCoreLinMeth.S = S;

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


%% REPLICATING BIGMACS MLN - Using WA2022 radiocarbon cores that were NOT in Lin2014 (and are >1000m and <40deg N or S latitudinally, and have >4 radiocarbon dates)

%Load the results again (just to ensure correct start point)
load("Results/dataT_PFnoLin_R200M20_Dec17.mat")

NoLinCoreRep.dataT = dataT;
NoLinCoreRep.S = S;

%Find number of cores
sizeCores = numel(NoLinCoreRep.dataT.lats);
numCores = sizeCores;

%Find high SR subset
highSRLog = NoLinCoreRep.dataT.meanSR > 8;

%Fit MLN distribution to my Bchron Mode data
disp("LinMethod (mode) data vs best fit mln")
[NoLinCoreRep.mixLogBMode, NoLinCoreRep.BModeHist,~,~,~,NoLinCoreRep.gmfitBmode, NoLinCoreRep.ncBmode, NoLinCoreRep.h, NoLinCoreRep.p,NoLinCoreRep.chiStat] = plotSRandResHistograms(dataT.bchronMode, x, highSRLog, 3, 1, 2, 0, "", 1, fitS);
if fitS.dispChi2
gcf; title("chi2gof: NoLinCoreBmode Data vs Best Fit MLN")
end

%Plot new histogram and MLN alongside BIGMACS MLN
figure;
hold on
yyaxis("left")
histogram(NoLinCoreRep.BModeHist, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts')
ylabel("Counts")
yyaxis("right")
plot(x, lognorm_BIGMACS.Var2, '-k', 'DisplayName', "BIGMACS", 'LineWidth', 1)
plot(x,NoLinCoreRep.mixLogBMode(:,2), '-r', 'DisplayName', "Replicated MLN", 'LineWidth', 1)
xlim([0 6])
xlabel("nSR")
ylabel("PDF")
legend()
title("BIGMACS replicated with independent cores")
set(gcf, 'Position', [20, 400, 500, 200])

%See how my distribution performs with chi2gof on BIGMACS gmfit
disp("LinMethod data (mode) vs BIGMACS MLN")
[NoLinCoreRep.hBM, NoLinCoreRep.pBM, NoLinCoreRep.chiStatBM] = chi2gof_vsMLN(gmfitBM, log(NoLinCoreRep.BModeHist), NoLinCoreRep.ncBmode, fitS);
if fitS.dispChi2
gcf; title("chi2gof: NoLinCoreBmode Data vs BIGMACS fit")
end

%Calculate TM of my Bchron Mode data
S.weighting = "age";
[~,~,NoLinCoreRep.TMnewBchNoWeight, NoLinCoreRep.TMnewBchWeightD] = TMcalculation(dataT.bchronMode, true(sizeCores), S);

%% REPLICATING BIGMACS MLN - Using WA2022 radiocarbon cores that were NOT in Lin2014 (and are >1000m and <40deg N or S latitudinally, and have >4 radiocarbon dates)

%Load the results again (just to ensure correct start point)
load("Results/dataT_PFnoLin_R200M20_Dec17.mat")

%Remove core(s)
logMD06_3075 = ismember(dataT.cores, 'MD06-3075');
NoLinCoreRep1.dataT = dataT(~logMD06_3075, :);
NoLinCoreRep1.S = S;

%Find number of cores
sizeCores = numel(NoLinCoreRep1.dataT.lats);
numCores = sizeCores;

%Fit MLN distribution to my Bchron Mode data
disp("LinMethod (mode) data vs best fit mln")
[NoLinCoreRep1.mixLogBMode, NoLinCoreRep1.BModeHist,~,~,~,NoLinCoreRep1.gmfitBmode, NoLinCoreRep1.ncBmode, NoLinCoreRep1.h, NoLinCoreRep1.p,NoLinCoreRep1.chiStat] = plotSRandResHistograms(NoLinCoreRep1.dataT.bchronMode, x, true(sizeCores), 3, 1, 2, 0, "", 1, fitS);
if fitS.dispChi2
gcf; title("chi2gof: NoLinCoreBmode Data vs Best Fit MLN")
end

%Plot new histogram and MLN alongside BIGMACS MLN
figure;
hold on
yyaxis("left")
histogram(NoLinCoreRep1.BModeHist, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts')
ylabel("Counts")
yyaxis("right")
plot(x, lognorm_BIGMACS.Var2, '-k', 'DisplayName', "BIGMACS", 'LineWidth', 1)
plot(x,NoLinCoreRep1.mixLogBMode(:,2), '-r', 'DisplayName', "Replicated MLN", 'LineWidth', 1)
xlim([0 6])
xlabel("nSR")
ylabel("PDF")
legend()
title("BIGMACS replicated with independent cores")
set(gcf, 'Position', [20, 400, 500, 200])

%See how my distribution performs with chi2gof on BIGMACS gmfit
disp("LinMethod data (mode) vs BIGMACS MLN")
[NoLinCoreRep1.hBM, NoLinCoreRep1.pBM, NoLinCoreRep1.chiStatBM] = chi2gof_vsMLN(gmfitBM, log(NoLinCoreRep1.BModeHist), NoLinCoreRep1.ncBmode, fitS);
if fitS.dispChi2
gcf; title("chi2gof: NoLinCoreBmode Data vs BIGMACS fit")
end

%The distribution here doesn't pass the chi2gof test for its own best-fit
%distribution. One core that contributes to this is probably MD06-3075,
%which has an nSR of 5.9 for 228 cm. Authors believe this is due to
%anomalous sediment oversampling during the coring process, as opposed to
%being due to a greater sedimentation rate. What happens without this core?

% %Calculate TM of my Bchron Mode data
% S.weighting = "age";
% [~,~,NoLinCoreRep1.TMnewBchNoWeight, NoLinCoreRep1.TMnewBchWeightD] = TMcalculation(dataT.bchronMode, true(sizeCores), S);
% 

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
    [h1R(i), p1R(i), chiStat1Run] = chi2gof_vsMLN(LinCoreLinMeth.gmfitBmode, log(LinCoreLinMeth.MLN1R.outputS.weightedC{i}), LinCoreLinMeth.MLN1R.outputS.numCpairs(i), fitS);
    chiStat1RunT = [chiStat1RunT; struct2table(chiStat1Run)]; %#ok<AGROW>
end

chiStat1RunT = addvars(chiStat1RunT, h1R, p1R, 'Before', "chi2stat");

LinCoreLinMeth.MLN1R.chiStat1RunT = chiStat1RunT;

% Test the chi2gof of each fit to BIGMACS
    h1R = NaN(numruns,1);
    p1R = NaN(numruns,1);
    fitS.dispChi2 = false;
    chiStat1RunT_BM = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);

for i = 1:numruns
    [h1R(i), p1R(i), chiStat1Run_BM] = chi2gof_vsMLN(gmfitBM, log(LinCoreLinMeth.MLN1R.outputS.weightedC{i}), LinCoreLinMeth.MLN1R.outputS.numCpairs(i), fitS);
    chiStat1RunT_BM = [chiStat1RunT_BM; struct2table(chiStat1Run_BM)]; %#ok<AGROW>
end

chiStat1RunT_BM = addvars(chiStat1RunT_BM, h1R, p1R, 'Before', "chi2stat");

LinCoreLinMeth.MLN1R.chiStat1RunT_BM = chiStat1RunT_BM;

%% Bring in results using updated dataset but still Marine09 and no reservoir uncertainty

load("Results/dataT_LinandPF_LinMethod_Dec10.mat")

NewCoreLinMeth.dataT = dataT;
NewCoreLinMeth.S = S;

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
    h1R = NaN(numruns,1);
    p1R = NaN(numruns,1);
    fitS.dispChi2 = false;
    chiStat1RunT = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);

for i = 1:numruns
    [h1R(i), p1R(i), chiStat1Run] = chi2gof_vsMLN(NewCoreLinMeth.gmfitBmode, log(NewCoreLinMeth.MLN1R.outputS.weightedC{i}), NewCoreLinMeth.MLN1R.outputS.numCpairs(i), fitS);
    chiStat1RunT = [chiStat1RunT; struct2table(chiStat1Run)]; %#ok<AGROW>
end

chiStat1RunT = addvars(chiStat1RunT, h1R, p1R, 'Before', "chi2stat");

NewCoreLinMeth.MLN1R.chiStat1RunT = chiStat1RunT;

% Test the chi2gof of each fit to BIGMACS
    h1R = NaN(numruns,1);
    p1R = NaN(numruns,1);
    fitS.dispChi2 = false;
    chiStat1RunT_BM = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);

for i = 1:numruns
    [h1R(i), p1R(i), chiStat1Run_BM] = chi2gof_vsMLN(gmfitBM, log(NewCoreLinMeth.MLN1R.outputS.weightedC{i}), NewCoreLinMeth.MLN1R.outputS.numCpairs(i), fitS);
    chiStat1RunT_BM = [chiStat1RunT_BM; struct2table(chiStat1Run)]; %#ok<AGROW>
end

chiStat1RunT_BM = addvars(chiStat1RunT_BM, h1R, p1R, 'Before', "chi2stat");

NewCoreLinMeth.MLN1R.chiStat1RunT_BM = chiStat1RunT_BM;

%% Bring in data with new dataset but with Marine 20 and 200 year R uncertainty

load("Results/dataT_LinandPF_R200M20_Dec12.mat")

NewCoreNewMeth.dataT = dataT;
NewCoreNewMeth.S = S;

sizeCores = numel(dataT.lats);
numCores = sizeCores;
fitS.dispChi2 = true;

%Fit MLN distribution to my Bchron Mode data
disp("LinMethod (mode) data vs best fit mln")
fitS.mlnReps = 5;
[NewCoreNewMeth.mixLogBMode, NewCoreNewMeth.BModeHist,~,~,~,NewCoreNewMeth.gmfitBmode, NewCoreNewMeth.ncBmode, NewCoreNewMeth.h, NewCoreNewMeth.p, NewCoreNewMeth.chiStat] = plotSRandResHistograms(dataT.bchronMode, x, true(sizeCores), 3, 1, 2, 0, "", 1, fitS);
if fitS.dispChi2
gcf; title("chi2gof: NewCoreBmode Data vs Best Fit MLN")
end
%See how my distributions perform with chi2gof on BIGMACS gmfit
disp("LinMethod data (mode) vs BIGMACS MLN")
[NewCoreNewMeth.hBM, NewCoreNewMeth.pBM, NewCoreNewMeth.chiStatBM] = chi2gof_vsMLN(gmfitBM, log(NewCoreNewMeth.BModeHist), NewCoreNewMeth.ncBmode, fitS);
if fitS.dispChi2
gcf; title("chi2gof: NewCoreBmode Data vs BIGMACS fit")
end
%Calculate TM of my Bchron Mode data
S.weighting = "age";
[~,~,NewCoreNewMeth.TMnewBchNoWeight, NewCoreNewMeth.TMnewBchWeightD] = TMcalculation(dataT.bchronMode, true(sizeCores), S);

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

%% Try with individual runs
numruns = 400;
fitS.dispChi2 = false;
[NewCoreNewMeth.MLN1R.pdfs, NewCoreNewMeth.MLN1R.c95up,...
    NewCoreNewMeth.MLN1R.c95down, NewCoreNewMeth.MLN1R.mus,...
    NewCoreNewMeth.MLN1R.sigmas, NewCoreNewMeth.MLN1R.outputS]...
    = SingleRunLogNorms(NewCoreNewMeth.dataT.bchronProb,...
    true(sizeCores, 1), numruns, x, 2, 3, 4, 0, fitS);

%%
% Test the chi2gof of each fit to BIGMACS
    h1R = NaN(numruns,1);
    p1R = NaN(numruns,1);
    fitS.dispChi2 = false;
    chiStat1RunT = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);

for i = 1:numruns
    [h1R(i), p1R(i), chiStat1Run] = chi2gof_vsMLN(NewCoreNewMeth.gmfitBmode, log(NewCoreNewMeth.MLN1R.outputS.weightedC{i}), NewCoreNewMeth.MLN1R.outputS.numCpairs(i), fitS);
    chiStat1RunT = [chiStat1RunT; struct2table(chiStat1Run)]; %#ok<AGROW>
end

chiStat1RunT = addvars(chiStat1RunT, h1R, p1R, 'Before', "chi2stat");

NewCoreNewMeth.MLN1R.chiStat1RunT = chiStat1RunT;

%%
% Test the chi2gof of each fit to BIGMACS
    h1R = NaN(numruns,1);
    p1R = NaN(numruns,1);
    fitS.dispChi2 = false;
    chiStat1RunT_BM = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);

for i = 1:numruns
    [h1R(i), p1R(i), chiStat1Run] = chi2gof_vsMLN(gmfitBM, log(NewCoreNewMeth.MLN1R.outputS.weightedC{i}), NewCoreNewMeth.MLN1R.outputS.numCpairs(i), fitS);
    chiStat1RunT_BM = [chiStat1RunT_BM; struct2table(chiStat1Run)]; %#ok<AGROW>
end

chiStat1RunT_BM = addvars(chiStat1RunT_BM, h1R, p1R, 'Before', "chi2stat");

NewCoreNewMeth.MLN1R.chiStat1RunT_BM = chiStat1RunT_BM;

%% Create table to show chi2Stats for each pairing
numPDFs = 4;
numDatas = 4;
chi2Table = table('Size', [numPDFs, numDatas], 'VariableTypes', ["double", "double", "double", "double"],'VariableNames',["BMhist", "LinCoresLinMethBMode", "NewCoresLinMethBMode", "NewCoresNewMethBMode"], 'RowNames',["BMpdf", "LCLMBModeMLN", "NCLMBModeMLN", "NCNMBModeMLN"]);
chi2Table.BMhist(1) = chiStatBMvBM.chi2stat;
chi2Table.LinCoresLinMethBMode(1) = LinCoreLinMeth.chiStatBM.chi2stat;
chi2Table.LinCoresLinMethBMode(2) = LinCoreLinMeth.chiStat.chi2stat;
chi2Table.NewCoresLinMethBMode(1) = NewCoreLinMeth.chiStatBM.chi2stat;
chi2Table.NewCoresLinMethBMode(3) = NewCoreLinMeth.chiStat.chi2stat;
chi2Table.NewCoresNewMethBMode(1) = NewCoreNewMeth.chiStatBM.chi2stat;
chi2Table.NewCoresNewMethBMode(4) = NewCoreNewMeth.chiStat.chi2stat;

figure;
hold on
plot(x, MLN_BIGMACS(:,2), '-k', 'DisplayName', "BIGMACS", 'LineWidth', 1)
plot(x, LinCoreLinMeth.mixLogBMode(:,2), '-b', 'DisplayName', "LinCores,Marine09", 'LineWidth', 1)
plot(x, NewCoreLinMeth.mixLogBMode(:,2), '-g', 'DisplayName', "NewCores,Marine09", 'LineWidth', 1)
plot(x, NewCoreNewMeth.mixLogBMode(:,2), '-r', 'DisplayName', "NewCores,Marine20", 'LineWidth', 1)
legend()
xlim([0 6])
xlabel("nSR")
ylabel("pdf")
title("Bchronology Mode Results")

%% Plot mix log normal from individual random sampling runs (note, random sampling of Bchronology at given depths - i.e. Lin Meth

commonYLim = [0 max([LinCoreLinMeth.MLN1R.c95up, NewCoreLinMeth.MLN1R.c95up, NewCoreNewMeth.MLN1R.c95up])];

figure;
hold on
plot(x, LinCoreLinMeth.MLN1R.pdfs, 'Color', [0,0,0,0.1], 'HandleVisibility', 'off')
plot(x, LinCoreLinMeth.mixLogBMode(:,2), '-r', 'DisplayName', "LinCores,Marine09", 'LineWidth', 1)
plot(NaN, NaN, 'DisplayName', "Random Runs that reject H0 = " + num2str(mean(LinCoreLinMeth.MLN1R.chiStat1RunT.h1R).*100, 3) + "%")
plot(x, MLN_BIGMACS(:,2), '-k', 'DisplayName', "BIGMACS", 'LineWidth', 1)
plot(NaN, NaN, 'DisplayName', "Random Runs that reject H0 = " + num2str(mean(LinCoreLinMeth.MLN1R.chiStat1RunT_BM.h1R).*100, 3) + "%")
xlim([0 6])
ylim(commonYLim)
legend()
title("LinCores; Marine09; R = 0±0")

figure;
hold on
plot(x, NewCoreLinMeth.MLN1R.pdfs, 'Color', [0,0,0,0.1], 'HandleVisibility', 'off')
plot(x, NewCoreLinMeth.mixLogBMode(:,2), '-r', 'DisplayName', "NewCores,Marine09", 'LineWidth', 1)
plot(NaN, NaN, 'DisplayName', "Random Runs that reject H0 = " + num2str(mean(NewCoreLinMeth.MLN1R.chiStat1RunT.h1R).*100, 3) + "%")
plot(x, MLN_BIGMACS(:,2), '-k', 'DisplayName', "BIGMACS", 'LineWidth', 1)
plot(NaN, NaN, 'DisplayName', "Random Runs that reject H0 = " + num2str(mean(NewCoreLinMeth.MLN1R.chiStat1RunT_BM.h1R).*100, 3) + "%")
xlim([0 6])
xlim([0 6])
ylim(commonYLim)
legend()
title("NewCores; Marine09; R = 0±0")

figure;
hold on
plot(x, NewCoreNewMeth.MLN1R.pdfs, 'Color', [0,0,0,0.1], 'HandleVisibility', 'off')
plot(x, NewCoreNewMeth.mixLogBMode(:,2), '-r', 'DisplayName', "BchronMode: NewCores,Marine20", 'LineWidth', 1)
plot(NaN, NaN,'LineStyle', "none", 'DisplayName', "Random Runs that reject H0 = " + num2str(mean(NewCoreNewMeth.MLN1R.chiStat1RunT.h1R).*100, 3) + "%")
plot(x, MLN_BIGMACS(:,2), '-k', 'DisplayName', "BchronMode: BIGMACS", 'LineWidth', 1)
plot(NaN, NaN, 'LineStyle', "none", 'DisplayName', "Random Runs that reject H0 = " + num2str(mean(NewCoreNewMeth.MLN1R.chiStat1RunT_BM.h1R).*100, 3) + "%")
xlim([0 6])
%ylim(commonYLim)
legend()
xlabel("nSR")
ylabel("PDF")
title("NewCores; Marine20; R = 0±200")