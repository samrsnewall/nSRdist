function[SR_MixLogNorm, histData, agediffsBinCounts, logSR_MixNorm, logSRbinCounts, gmfit, numSRcalcs, totalSedLength, totalSedAge] = plotSRandResHistograms(nSRcounts, x, coreSubsetLogical, weightRepDP, weightRepInflator, components, regularizationValue, subsetName, plotQ, fitS)
%%% This function takes some normalised sedimentation rate data in cell
%%% format, combines the counts into a single array, applies the weighting
%%% and then fits a mixture log normal to the result. It will also plot the
%%% histogram of the logarithm data with the mixture normal on top, show a
%%% QQ plot of the logarithm data, show the histogram of the data with the
%%% log mix normal on top, and plot a histogram of the resolution of the
%%% data

%% Combine all counts into one array
%set up arrays to be concatenated into
nSRcountsArraywNaN = countsCell2Array(nSRcounts, coreSubsetLogical);

%Remove NaNs that are used to separate cores and runs
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

data = makeWeightedReplicates(nSR, weightingsArray, weightRepDP, weightRepInflator);

numSRcalcs = size(nSRcountsArray, 2);
totalSedLength = sum(depthDiffs, "omitnan");
totalSedAge = sum(ageDiffs, "omitnan");

%Calculate the MixLogNorm from these counts
[a,b] = size(data);
if a>b
    data = data';
end
[SR_MixLogNorm, logSR_MixNorm, gmfit] = fitMixLogNorm(data, x, components, fitS.mlnReps);

fitInfo.nll = gmfit.NegativeLogLikelihood;
fitInfo.BIC = gmfit.BIC;

dataLog = log(data);

% %% Perform chi2gof on gmfit
% %Create cdf of distribution
% %Get the parameters of the mixed gaussian in log space
% mu1 = gmfit.mu(1); sigma1 = sqrt(gmfit.Sigma(1));
% mu2 = gmfit.mu(2); sigma2 = sqrt(gmfit.Sigma(2));
% w1  = gmfit.ComponentProportion(1);
% w2  = gmfit.ComponentProportion(2);
% 
% %Create the cdf function handle
% mlncdf = @(t) w1 * normcdf(t, mu1, sigma1) + w2 * normcdf(t, mu2, sigma2);

%% count how many estimates of nSR
numbernSRcounts = length(nSR);

%% Calculate mean and variance of data
nSRdataAsCounts = data;
nSRcutoff = nSRdataAsCounts > 10;
dataMeanLinear = geomean(nSRdataAsCounts);
dataVarLinear = var(nSRdataAsCounts(~nSRcutoff));
dataMeanLog = mean(dataLog);
dataVarLog = var(dataLog);

%% Calculate Histogram Bin Counts
%Define bins edges for histogram
SRbinEdges = 0:0.1:10;
logSRbinEdges = -4:0.05:4;
agediffsBinEdges = 0:100:10000;

%Calculate bin counts (using weighting for SR)
%SRbinCounts = makeWeightedBinCounts(nSRcountsArray(1,:), nSRcountsArray(2,:), SRbinEdges);
histData = data;
logSRbinCounts = makeWeightedBinCounts(log(nSR), depthDiffs, logSRbinEdges);
agediffsBinCounts = makeWeightedBinCounts(ageDiffs, depthDiffs, agediffsBinEdges);

%% Compare histogram of "weighted data" to weighted histogram of data
% figure;
% subplot(1,2,1)
% histogram("BinCounts", SRbinCounts, "BinEdges", SRbinEdges)
% xlim([0 6])
% title("Weighted by summing weights")
% subplot(1,2,2)
% histogram(data, "BinEdges", SRbinEdges)
% xlim([0 6])
% title("Weighted by replicating data")
%
% %% compare histogram of weighted log data to weighted histogram of log data
% figure;
% subplot(1,2,1)
% histogram("BinCounts", logSRbinCounts, "BinEdges", logSRbinEdges)
% xlim([-4 4])
% title("Weighted by summing weights")
% subplot(1,2,2)
% histogram(dataLog, "BinEdges", logSRbinEdges)
% xlim([-4 4])
% title("Weighted by replicating data")
%
% %% Get true weighted histograms of nSR and lognSR
% figure;
% subplot(2,2,1)
% histogram("BinCounts", logSRbinCounts, "BinEdges", logSRbinEdges, "Normalization", "pdf")
% xlim([-4 4])
% title("Weighted by summing weights")
% subplot(2,2,2)
% histogram("BinCounts", SRbinCounts, "BinEdges", SRbinEdges, "Normalization", "pdf")
% xlim([0 6])
% sgtitle("True histograms")
% title("Weighted by summing weights")
% subplot(2,2,3)
% histogram(dataLog, "BinEdges", logSRbinEdges, "Normalization", "pdf")
% xlim([-4 4])
% title("Weighted by replicating data")
% subplot(2,2,4)
% histogram(data, "BinEdges", SRbinEdges, "Normalization", "pdf")
% xlim([0 6])
% title("Weighted by replicating data")


%% Plot to see how data compare to estimated distributions
if plotQ == 1
    figure
    subplot(2,2,1)
    hold on
    histogram(dataLog, "Normalization", "pdf")
    plot(logSR_MixNorm(:,1), logSR_MixNorm(:,2))
    xlabel("Log nSR")
    title("Data Mean = " + num2str(dataMeanLog) + "; Data Var = " + num2str(dataVarLog))
    subplot(2,2,2)
    qqplot(dataLog)
    xlabel("Log nSR Data")
    subplot(2,2,3)
    hold on
    histogram(data, "BinEdges", SRbinEdges, "Normalization", "pdf")
    plot(SR_MixLogNorm(:,1), SR_MixLogNorm(:,2))
    xlabel("nSR")
    xlim([0 6])
    title("Data Mean = " + num2str(dataMeanLinear) + "; Data Var = " + num2str(dataVarLinear))
    subplot(2,2,4)
    histogram("BinCounts", agediffsBinCounts, "BinEdges", agediffsBinEdges)
    xlabel("Age Diff (yrs)")
    ylabel("Counts")
    title(num2str(numbernSRcounts) + " nSR estimates")
    sgtitle(subsetName)
end
%
% figure
% subplot(1,2,1)
% hold on
% histogram("BinCounts", SRbinCounts, "BinEdges", SRbinEdges, 'FaceColor', 'k', 'FaceAlpha', 0.5, "Normalization","pdf")
% plot(SR_MixLogNorm(:,1), SR_MixLogNorm(:,2), 'Color', 'k', 'LineWidth', 2)
% xlabel("nSR")
% xlim([0 6])
% title("Data Mean = " + num2str(dataMeanLinear) + "; Data Var = " + num2str(dataVarLinear))
% subplot(1,2,2)
% histogram("BinCounts", agediffsBinCounts, "BinEdges", agediffsBinEdges,'FaceColor', 'k', 'FaceAlpha', 0.5)
% xlabel("Age Diff (yrs)")
% ylabel("Length of sediment (cm)")
% sgtitle(subsetName)


end