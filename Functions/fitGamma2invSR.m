function[nSR_invGamma, nSR_invGammaProb, invxGamProb, alpha, invSRbincounts_weighted] = fitGamma2invSR(nSRcounts, coreSubsetLogical,weightDP, weightInflator, x, fitS)
% This function takes nSRcounts data, inverts it so that it is accumulation
% rate data (as discussed in Blaauw et al., 2011) and fits a gamma
% distribution to the data. It then inverts that gamma distribution so that
% it is over nSR again.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine all counts into one array
nSRcountsArraywNaN = countsCell2Array(nSRcounts, coreSubsetLogical);

% Clean up data
removeIndex = isnan(nSRcountsArraywNaN(1,:)); %Get rid of nans
removeIndex2 = nSRcountsArraywNaN(1,:) == 0; %Get rid of zeros
nSRcountsArray = nSRcountsArraywNaN(:, ~(removeIndex | removeIndex2)); % remove nans and zeros

%Remove information from age pairs not within range from
%fitS.Lin2014AgeFilter
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

% Invert the data so that it's in inverse sed rate
invnSRcounts = 1./nSRcountsArray(1,:);

%Get weighted replicates
invnSR_WR = makeWeightedReplicates(invnSRcounts, weightingsArray, weightDP, weightInflator);

%Create weighted bin counts
invSRbincounts_weighted = makeWeightedBinCounts(invnSRcounts, weightingsArray, fitS.invXbinEdges);

% Set the range of x values of interest 
% (use x values currently used in BIGMACS)
invx = sort(1./x);

% Estimate the gamma fit parameters
phat = gamfit(invnSR_WR);
alpha = phat(1);
beta = phat(2);

% Create gamma on invx values
invxGamProb = gampdf(invx, alpha, beta);

% Convert back to nSR
f_inv = @(x) 1./x;
[nSR_invGamma, nSR_invGammaProb] = px_to_pfx(invx, invxGamProb, f_inv);
end