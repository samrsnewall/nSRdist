%% Bring in results using updated dataset but still Marine09 and no reservoir uncertainty
%Add important paths
addpath('../Functions')

%Load BIGMACS files
lognorm_BIGMACS = readtable("../BIGMACSdata/lognormal.txt");                  % (use x values currently used in BIGMACS)
x = lognorm_BIGMACS.Var1';
BIGMACShist = readmatrix("../BIGMACSdata/Lin2014_sedrateratio_cm_wo_NaN.txt");
BIGMACShist = BIGMACShist(:,4);
TM_BIGMACS = readmatrix("../BIGMACSdata/transition_parameter.txt");

%Load results 
load("../Results/dataT_LinandPF_LinMethod_Dec10.mat")
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
