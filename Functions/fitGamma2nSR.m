function[nSR_Gamma, nSR_GammaProb, alpha, nSRbincounts_weighted] = fitGamma2nSR(nSRcounts, coreSubsetLogical, weightDP, weightInflator, x, fitS)
% This function takes nSRcounts data and fits a gamma
% distribution to the data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine all counts into one array
nSRcountsArraywNaN = countsCell2Array(nSRcounts, coreSubsetLogical);

% Clean up data
removeIndex = isnan(nSRcountsArraywNaN(1,:)); %Get rid of nans
removeIndex2 = nSRcountsArraywNaN(1,:) == 0; %Get rid of zeros
nSRcountsArray = nSRcountsArraywNaN(:, ~(removeIndex | removeIndex2)); % remove nans and zeros

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

%Apply weighting
if fitS.weighting == "none"
    weightingsArray = ones(1,length(nSRcountsArray(2,:)));
elseif fitS.weighting == "depth"
    weightingsArray = nSRcountsArray(3,:);
elseif fitS.weighting == "age"
    weightingsArray = nSRcountsArray(4,:);
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
alpha = phat(1);
beta = phat(2);

% Create gamma
nSR_GammaProb = gampdf(x, alpha, beta);
end