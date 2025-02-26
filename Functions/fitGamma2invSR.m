function[nSR_invGamma, nSR_invGammaProb, invxGamProb, alpha, invSRbincounts_weighted] = fitGamma2invSR(nSRcounts, coreSubsetLogical, fitS)
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

%Remove information from age pairs not within 0.5-4.5 kyr
if fitS.Lin2014AgeFiltering
    if sum(nSRcountsArray(4,:) < 0) ~= 0
        warning("There are negative sed rates being filtered out!")
    end
    L2014Log = nSRcountsArray(4,:) < 4500 & nSRcountsArray(4,:) > 500;
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

% Invert the data so that it's in inverse sed rate
invnSRcounts = 1./nSRcountsArray(1,:);

%Get weighted replicates
invnSR_WR = makeWeightedReplicates(invnSRcounts, weightingsArray, 3, 10);

%Create weighted bin counts
invSRbincounts_weighted = makeWeightedBinCounts(invnSRcounts, weightingsArray, fitS.invXbinEdges);

% Set the range of x values of interest 
% (use x values currently used in BIGMACS)
lognorm_BIGMACS = readtable("../lognormal_BIGMACS.txt");
x = lognorm_BIGMACS.Var1';
invx = sort(1./x);

% Estimate the gamma fit parameters
phat = gamfit(invnSR_WR);
alpha = phat(1);
beta = phat(2);

% Create gamma on invx values
invxGamProb = gampdf(invx, alpha, beta);
% invxGamcdf = @(t) gamcdf(t, alpha, beta);
% 
% %Set up certain important parameters
% %desiredSum = length(dataLog);
% desiredSum = numSRcalcs;
% binN = fitS.chi2binN;

% %Perform chi2gof
% numParams = 2;
% [h,p,chiStat] = chi2gof_vsfunction(invnSR_WR, invxGamcdf, numParams, desiredSum, binN, fitS);
% gcf;
% title("chi2gof of Data vs Best Fit Gamma")

% % Plot inverse data
% binEdges = [0:0.1:15];
% 
% figure;
% hold on
% histogram(invnSR_WR, "BinEdges", binEdges, 'Normalization', "pdf");
% plot(invx, invxGamProb, 'LineWidth', 2)
% xlim([0 5])
% xlabel("Inverse normalised SR")
% title("\alpha = " + num2str(phat(1)))

% Convert back to nSR
[nSR_invGamma, nSR_invGammaProb] = gammaAccRate2nSR(phat(1), phat(2));
end