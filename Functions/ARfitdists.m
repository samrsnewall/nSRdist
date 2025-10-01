function[dStru] = ARfitdists(dataTCol, x, chooseLog, weightDP, weightInflator, countDivisor, fitS)
%This function fits the 4 distributions to all samples in a dataT column
%(used for BMode and can be used to fit dists to all BSamp samples instead
%of fitting to each of 1000 BSamp runs)

%Fit Mix Log Normal
[dStru.mixLog, dStru.weightedC,dStru.agediffs,~,~,...
    dStru.MLN.nSR.gmfit, dStru.numCpairs, dStru.sedLength, ...
    dStru.sedTimeSpan]...
    = plotSRandResHistograms(dataTCol,...
    x, chooseLog, weightDP, weightInflator, 2, 0, "", 0, fitS);

dStru.numCpairs = dStru.numCpairs./countDivisor;
dStru.MLN.nSR.x = dStru.mixLog(:,1);
dStru.MLN.nSR.px = dStru.mixLog(:,2);
[dStru.MLN.nSR.mu, dStru.MLN.nSR.var] = muVarPDFVec(dStru.MLN.nSR);
dStru.MLN.nSR.fitInfo.nll = dStru.MLN.nSR.gmfit.NegativeLogLikelihood;
dStru.MLN.nSR.fitInfo.BIC = dStru.MLN.nSR.gmfit.BIC;

%Fit Inverse gamma distribution to data
[dStru.invGam.nSR.x, dStru.invGam.nSR.px, ~,~,~,dStru.invGam.nSR.fitInfo] = fitGamma2invSR(dataTCol, chooseLog, weightDP, weightInflator, x, fitS);
[dStru.invGam.nSR.mu, dStru.invGam.nSR.var] = muVarPDFVec(dStru.invGam.nSR);

%Fit gamma distribution to data
[dStru.Gam.nSR.x, dStru.Gam.nSR.px, dStru.LN.nSR.fitInfo] = fitGamma2nSR(dataTCol, chooseLog, weightDP, weightInflator, x, fitS);
[dStru.Gam.nSR.mu, dStru.Gam.nSR.var] = muVarPDFVec(dStru.Gam.nSR);

%Fit LogNormal distribution to data
[dStru.LN.nSR.x, dStru.LN.nSR.px, dStru.LN.nSR.fitInfo] = fitLogNorm2nSR(dataTCol, chooseLog, weightDP, weightInflator, x, fitS);
[dStru.LN.nSR.mu, dStru.LN.nSR.var] = muVarPDFVec(dStru.LN.nSR);

%% Transform pdf Vecs to log space
f_log = @(x) log(x);
[dStru.MLN.lnSR.x, dStru.MLN.lnSR.px] = px_to_pfx(dStru.MLN.nSR.x, dStru.MLN.nSR.px, f_log);
[dStru.MLN.lnSR.mu, dStru.MLN.lnSR.var] = muVarPDFVec(dStru.MLN.lnSR);
dStru.MLN.lnSR.numParams = 5;
dStru.MLN.lnSR.pdfName = "2 Component Mix Log Normal";

[dStru.invGam.lnSR.x, dStru.invGam.lnSR.px] = px_to_pfx(dStru.invGam.nSR.x, dStru.invGam.nSR.px, f_log);
[dStru.invGam.lnSR.mu, dStru.invGam.lnSR.var] = muVarPDFVec(dStru.invGam.lnSR);
dStru.invGam.lnSR.numParams = 2;
dStru.invGam.lnSR.pdfName = "Inverse Gamma";

[dStru.Gam.lnSR.x, dStru.Gam.lnSR.px] = px_to_pfx(dStru.Gam.nSR.x, dStru.Gam.nSR.px, f_log);
[dStru.Gam.lnSR.mu, dStru.Gam.lnSR.var] = muVarPDFVec(dStru.Gam.lnSR);
dStru.Gam.lnSR.numParams = 2;
dStru.Gam.lnSR.pdfName = "Gamma";

[dStru.LN.lnSR.x, dStru.LN.lnSR.px] = px_to_pfx(dStru.LN.nSR.x, dStru.LN.nSR.px, f_log);
[dStru.LN.lnSR.mu, dStru.LN.lnSR.var] = muVarPDFVec(dStru.LN.lnSR);
dStru.LN.lnSR.numParams = 2;
dStru.LN.lnSR.pdfName = "LogNorm";

%Make a column cell-vector that holds all pdfs to test
pdfs = {dStru.LN.lnSR; dStru.MLN.lnSR; dStru.Gam.lnSR; dStru.invGam.lnSR};
%pdfs = {dStru.MLN.lnSR; dStru.invGam.lnSR};
fitS.dispChi2 = true;
% h = NaN(1,1);
% p = NaN(1,1);
[h, p, chiStat] = chi2_dataVStwopdfVECs(log(dStru.weightedC), dStru.numCpairs, 20, pdfs, fitS);

%store the chiStats in the data structures
for i = 1:size(pdfs,1)
    chiStat{i}.h = h(i);
    chiStat{i}.p = p(i);
end
dStru.LN.chiStats = chiStat{1};
dStru.MLN.chiStats = chiStat{2};
dStru.Gam.chiStats = chiStat{3};
dStru.invGam.chiStats = chiStat{4};
end