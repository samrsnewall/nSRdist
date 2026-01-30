function[nSR_Gamma, nSR_GammaProb, fitInfo, nSRbincounts_weighted] = fitGamma2nSR(nSRcounts, coreSubsetLogical, weightDP, weightInflator, x, fitS)
% This function takes nSRcounts data and fits a gamma
% distribution to the data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine all counts into one array
nSRcountsArraywNaN = countsCell2Array(nSRcounts, coreSubsetLogical);

% Clean up data
NaNLog = isnan(nSRcountsArraywNaN(1,:)); %Get rid of nans
ZerosLog = nSRcountsArraywNaN(1,:) == 0; %Get rid of zeros
nSRcountsArray = nSRcountsArraywNaN(:, ~(NaNLog | ZerosLog)); % remove nans and zeros

%Remove information from age pairs not within 0.5-4 kyr
if fitS.Lin2014AgeFiltering
    if sum(nSRcountsArray(4,:) < 0) ~= 0
        warning("There are negative sed rates being filtered out!")
    end
    L2014Log = nSRcountsArray(4,:) < max(fitS.Lin2014AgeFilter) & nSRcountsArray(4,:) > min(fitS.Lin2014AgeFilter);
    nSRcountsArray = nSRcountsArray(:,L2014Log);
else
    if sum(nSRcountsArray(4,:) < 0) ~= 0
        error("There are negative sed rates")
    end
end

%Convert the weighted nSR counts to a single dimension array of counts
%where the number of counts is representative of their weighting
nSR = nSRcountsArray(1,:)'; %nSR data
depthDiffs = nSRcountsArray(3,:);  %weightings
ageDiffs = nSRcountsArray(4,:);
if fitS.weighting == "none"
    data = nSR;
    weightingsArray = ones(1,length(nSRcountsArray(2,:)));
elseif fitS.weighting == "depth"
        weightingsArray = depthDiffs;
elseif fitS.weighting == "age"
        weightingsArray = ageDiffs;
end
numSRcalcs = size(nSRcountsArray, 2);
nSRcounts = nSRcountsArray(1,:);

%Get weighted replicates
nSR_WR = makeWeightedReplicates(nSRcounts, weightingsArray, weightDP, weightInflator);

%Create weighted bin counts
nSRbincounts_weighted = makeWeightedBinCounts(nSRcounts, weightingsArray, fitS.invXbinEdges);

nSR_Gamma = x;

% Estimate the gamma fit parameters
phat = gamfit(nSR_WR);
fitInfo.alpha = phat(1);
fitInfo.beta = phat(2);
fitInfo.nll = gamlike([fitInfo.alpha, fitInfo.beta], nSR_WR);
fitInfo.NumParams = 2;
fitInfo.BIC = 2*fitInfo.nll + fitInfo.NumParams*log(length(nSR_WR));

% Create gamma
nSR_GammaProb = gampdf(x, fitInfo.alpha, fitInfo.beta);
end