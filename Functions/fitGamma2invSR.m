function[nSR_invGamma, nSR_invGammaProb, invxGamProb, alpha, invSRbincounts_weighted, fitInfo] = fitGamma2invSR(nSRcounts, coreSubsetLogical,weightDP, weightInflator, x, fitS)
% This function takes nSRcounts data, inverts it so that it is accumulation
% rate data (as discussed in Blaauw et al., 2011) and fits a gamma
% distribution to the data. It then inverts that gamma distribution so that
% it is over nSR again.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Combine all counts into one array
%set up arrays to be concatenated into
nSRcountsArraywNaN = countsCell2Array(nSRcounts, coreSubsetLogical);

% Clean up data
NaNLog = isnan(nSRcountsArraywNaN(1,:)); %Get rid of nans
ZerosLog = nSRcountsArraywNaN(1,:) == 0; %Get rid of zeros
nSRcountsArray = nSRcountsArraywNaN(:, ~(NaNLog | ZerosLog)); % remove nans and zeros

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

%% Fit Mix Log Norm
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
fitInfo.nll = gamlike([alpha, beta], invnSR_WR);

% Create gamma on invx values
invxGamProb = gampdf(invx, alpha, beta);

% Convert back to nSR
f_inv = @(x) 1./x;
[nSR_invGamma, nSR_invGammaProb] = px_to_pfx(invx, invxGamProb, f_inv);
end