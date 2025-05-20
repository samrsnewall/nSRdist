function[nSR_LogNorm, nSR_LogNormProb, gmfit, nSRbincounts_weighted] = fitLogNorm2nSR(nSRcounts, coreSubsetLogical, fitS)
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
    weightingsArray = nSRcountsArray(2,:);
elseif fitS.weighting == "age"
    weightingsArray = nSRcountsArray(4,:);
end
numSRcalcs = size(nSRcountsArray, 2);
nSRcounts = nSRcountsArray(1,:);

%Get weighted replicates
nSR_WR = makeWeightedReplicates(nSRcounts, weightingsArray, 3, 10);

%Create weighted bin counts
nSRbincounts_weighted = makeWeightedBinCounts(nSRcounts, weightingsArray, fitS.invXbinEdges);

% Set the range of x values of interest 
% (use x values currently used in BIGMACS)
lognorm_BIGMACS = readtable("../lognormal_BIGMACS.txt");
x = lognorm_BIGMACS.Var1';

%Fit mix log norm with 1 component
[nSR_LogNormVec, logSR_NormVec, gmfit] = fitMixLogNorm(nSR_WR, x, 1, fitS.mlnReps);

nSR_LogNorm = nSR_LogNormVec(:,1);
nSR_LogNormProb = nSR_LogNormVec(:,2);

end