function[nSR_Gamma, nSR_GammaProb, invxGamProb, alpha, invSRbincounts_weighted] = fitGamma2invSR(nSRcounts, coreSubsetLogical, invSRbinedges)
% This function takes nSRcounts data, inverts it so that it is accumulation
% rate data (as discussed in Blaauw et al., 2011) and fits a gamma
% distribution to the data. It then inverts that gamma distribution so that
% it is over nSR again.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine all counts into one array
%set up arrays to be concatenated into
nSRcountsArrayRAW = ones(4,1);
%Concatenate all nSRcounts that are in the desired subset
for i = 1:length(nSRcounts)
    if coreSubsetLogical(i) == 1
    nSRcountsArrayRAW = cat(2, nSRcountsArrayRAW, nSRcounts{i});
    end
end
%Remove the ones that were used to set up arrays
nSRcountsArraywNaN = nSRcountsArrayRAW(:,2:end);

% Clean up data
removeIndex = isnan(nSRcountsArraywNaN(1,:)); %Get rid of nans
removeIndex2 = nSRcountsArraywNaN(1,:) == 0; %Get rid of zeros
nSRcountsArray = nSRcountsArraywNaN(1, ~(removeIndex | removeIndex2)); % remove nans and zeros
weightingsArray = nSRcountsArraywNaN(2,~(removeIndex | removeIndex2));

% Invert the data so that it's in inverse sed rate
invnSRcounts = 1./nSRcountsArray;

%Get weighted replicates
invnSR_WR = makeWeightedReplicates(invnSRcounts, weightingsArray, 3, 10);

%Create weighted bin counts
invSRbincounts_weighted = makeWeightedBinCounts(invnSRcounts, weightingsArray, invSRbinedges);

% Set the range of x values of interest 
% (use x values currently used in BIGMACS)
lognorm_BIGMACS = readtable("lognormal_BIGMACS.txt");
x = lognorm_BIGMACS.Var1';
invx = sort(1./x);


% Estimate the gamma fit parameters
phat = gamfit(invnSR_WR);
alpha = phat(1);
beta = phat(2);

% Create gamma on invx values
invxGamProb = gampdf(invx, alpha, beta);

% Plot inverse data
binEdges = [0:0.1:15];

figure;
hold on
histogram(invnSR_WR, "BinEdges", binEdges, 'Normalization', "pdf");
plot(invx, invxGamProb, 'LineWidth', 2)
xlim([0 5])
xlabel("Inverse normalised SR")
title(num2str(phat(1)))

% Convert back to nSR
[nSR_Gamma, nSR_GammaProb] = gammaAccRate2nSR(phat(1));

end