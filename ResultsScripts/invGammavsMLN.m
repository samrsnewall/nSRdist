%% Bring in data with new dataset but with Marine 20 and 200 year R uncertainty
%Add important paths
addpath('../Functions')

%Load BIGMACS files
lognorm_BIGMACS = readtable("../BIGMACSdata/lognormal.txt");                  % (use x values currently used in BIGMACS)
x = lognorm_BIGMACS.Var1';
BIGMACShist = readmatrix("../BIGMACSdata/Lin2014_sedrateratio_cm_wo_NaN.txt");
BIGMACShist = BIGMACShist(:,4);
TM_BIGMACS = readmatrix("../BIGMACSdata/transition_parameter.txt");

%Load results file
load("../Results/dataT_LinandPF_LinMethod_Dec10.mat")
NewCoreNewMeth.dataT = dataT;
NewCoreNewMeth.S = S;
sizeCores = numel(dataT.lats);
numCores = sizeCores;

%Set up fitS structure
fitS.dispChi2 = true;
fitS.Lin2014AgeFiltering = true;
fitS.mlnReps = 5;
fitS.weighting = "depth";
fitS.chi2binN = 10;
fitS.invXbinEdges = 0:0.1:15;
fitS.enforceBinSizeLimits = false;

% Want to fit some inverse gamma distributions and mix log normal
% distributions to the same data to see how they look, and compare their
% chi2stat values...

%Fit MLN distribution to my Bchron Mode data
disp("LinMethod (mode) data vs best fit mln")
[NewCoreNewMeth.mixLogBMode, NewCoreNewMeth.BModeHist,~,~,~,...
    NewCoreNewMeth.gmfitBmode, NewCoreNewMeth.ncBmode, NewCoreNewMeth.h,...
    NewCoreNewMeth.p, NewCoreNewMeth.chiStat]...
    = plotSRandResHistograms(dataT.bchronMode, x, true(sizeCores), 3, 1,...
    2, 0, "", 1, fitS);
if fitS.dispChi2
    gcf; title("chi2gof: NewCoreBmode Data vs Best Fit MLN")
end

%Fit Inverse gamma distribution to my Bchron Mode data
[nSR_Gamma, nSR_GammaProb] = fitGamma2invSR(NewCoreNewMeth.dataT.bchronMode, true(sizeCores), fitS);

%Plot 


