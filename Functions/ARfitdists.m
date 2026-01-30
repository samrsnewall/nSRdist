function[outS] = ARfitdists(dataTCol, x, chooseLog, weightRepDP, weightRepInflator, countDivisor, fitS)
%This function fits the 4 distributions to all samples in a dataT column
%(used for BMode and can be used to fit dists to all BSamp samples instead
%of fitting to each of 1000 BSamp runs)

%% Combine all counts into one array
%set up arrays to be concatenated into
nSRcountsArray = countsCell2Array(dataTCol, chooseLog);



%Get rid of negative sedimentation rates (ideally should not be there)
if sum(nSRcountsArray(1,:) < 0) ~= 0
    warning("There are negative sed rates!... being removed")
    nSRcountsArray = nSRcountsArray(:, nSRcountsArray(1,:) >= 0);
end

%Get rid of SRs of 0 (ideally shouldn't be there)
if sum(nSRcountsArray(1,:) == 0) ~= 0
    warning("There are SRs of 0!... being removed")
    nSRcountsArray = nSRcountsArray(:, nSRcountsArray(1,:) > 0);
end

%Apply filtering as done by Lin2014 if desired
if fitS.Lin2014AgeFiltering
    L2014Log = nSRcountsArray(4,:) < max(fitS.Lin2014AgeFilter) & nSRcountsArray(4,:) > min(fitS.Lin2014AgeFilter);
    nSRcountsArray = nSRcountsArray(:,L2014Log);
end

%% Apply weighting
%Convert the weighted nSR counts to a single dimension array of counts
%where the number of counts is representative of their weighting

%Remove NaNs that are used to separate cores and runs

inputData = nSRcountsArray(1,:); %nSR data
inputData = inputData(~isnan(inputData)); %nSR data
depthDiffs = nSRcountsArray(3,:);  %weightings
ageDiffs = nSRcountsArray(4,:);
if fitS.weighting == "none"
    data = inputData;
    %weightingsArray = ones(1,length(nSRcountsArray(2,:)));
elseif fitS.weighting == "depth"
        weightingsArray = depthDiffs;
        data = makeWeightedReplicates(inputData, weightingsArray, weightRepDP, weightRepInflator);
elseif fitS.weighting == "age"
        weightingsArray = ageDiffs;
        data = makeWeightedReplicates(inputData, weightingsArray, weightRepDP, weightRepInflator);
end

outS.numCpairs = size(nSRcountsArray, 2)./countDivisor;

%% Resample data if desired
if fitS.resampleData
    rng(3)
    alldata = data;
    data = randsample(alldata, floor(outS.numCpairs), false);
end

%% Calculate the MixLogNorm from these counts
% [a,b] = size(data);
% if a>b
%     data = data';
% end
[SR_MixLogNorm, ~, MLN.fitInfo] = fitMixLogNorm(data, x, 2, fitS.mlnReps);
MLN.nSR.x = SR_MixLogNorm(:,1);
MLN.nSR.px = SR_MixLogNorm(:,2);

%% Fit Lognormal to these counts
[SR_LogNorm, ~, LN.fitInfo] = fitMixLogNorm(data, x, 1, fitS.mlnReps);
LN.nSR.x = SR_LogNorm(:,1);
LN.nSR.px = SR_LogNorm(:,2);

%% Fit gamma to these counts
% Estimate the gamma fit parameters
[Gam.nSR.px, Gam.fitInfo] = fitGamma(data, x);
Gam.nSR.x = x';

%% Fit inverse gamma to these counts (i.e. fit gamma to inverse of counts)
[invGam.nSR.px, ~, invGam.fitInfo] = fitInvGamma(data, x);
invGam.nSR.x = x';


%% Save other useful information

outS.numCpairs = size(nSRcountsArray, 2)./countDivisor;
outS.weightedC = data;

agediffsBinEdges = 0:100:10000;
outS.agediffsWC= makeWeightedBinCounts(ageDiffs, depthDiffs, agediffsBinEdges);

outS.sedLength = sum(depthDiffs, "omitnan");
outS.sedTimeSpan = sum(ageDiffs, "omitnan");

%Calculate mean and variance of each pdf
[invGam.nSR.mu, invGam.nSR.var] = muVarPDFVec(invGam.nSR);
[Gam.nSR.mu, Gam.nSR.var] = muVarPDFVec(Gam.nSR);
[LN.nSR.mu, LN.nSR.var] = muVarPDFVec(LN.nSR);
[MLN.nSR.mu, MLN.nSR.var] = muVarPDFVec(MLN.nSR);


%%


%% Transform pdf Vecs to log space and set up fo
f_log = @(x) log(x);
[MLN.lnSR.x, MLN.lnSR.px] = px_to_pfx(MLN.nSR.x, MLN.nSR.px, f_log);
[MLN.lnSR.mu, MLN.lnSR.var] = muVarPDFVec(MLN.lnSR);
MLN.lnSR.numParams = 5;
MLN.lnSR.pdfName = "2 Component Mix Log Normal";

[invGam.lnSR.x, invGam.lnSR.px] = px_to_pfx(invGam.nSR.x, invGam.nSR.px, f_log);
[invGam.lnSR.mu, invGam.lnSR.var] = muVarPDFVec(invGam.lnSR);
invGam.lnSR.numParams = 2;
invGam.lnSR.pdfName = "Inverse Gamma";

[Gam.lnSR.x, Gam.lnSR.px] = px_to_pfx(Gam.nSR.x, Gam.nSR.px, f_log);
[Gam.lnSR.mu, Gam.lnSR.var] = muVarPDFVec(Gam.lnSR);
Gam.lnSR.numParams = 2;
Gam.lnSR.pdfName = "Gamma";

[LN.lnSR.x, LN.lnSR.px] = px_to_pfx(LN.nSR.x, LN.nSR.px, f_log);
[LN.lnSR.mu, LN.lnSR.var] = muVarPDFVec(LN.lnSR);
LN.lnSR.numParams = 2;
LN.lnSR.pdfName = "LogNorm";

%Make a column cell-vector that holds all pdfs to test
pdfs = {LN.lnSR; MLN.lnSR; Gam.lnSR; invGam.lnSR};
%pdfs = {MLN.lnSR; invGam.lnSR};
fitS.dispChi2 = false;
% h = NaN(1,1);
% p = NaN(1,1);
[h, p, chiStat] = chi2_dataVStwopdfVECs(log(data), outS.numCpairs, 20, pdfs, fitS);

%store the chiStats in the data structures
for i = 1:size(pdfs,1)
    chiStat{i}.h = h(i);
    chiStat{i}.p = p(i);
end
LN.chiStats = chiStat{1};
MLN.chiStats = chiStat{2};
Gam.chiStats = chiStat{3};
invGam.chiStats = chiStat{4};

%% Put structures into output structure
outS.LN = LN;
outS.MLN = MLN;
outS.Gam = Gam;
outS.invGam = invGam;


end