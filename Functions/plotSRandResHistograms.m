function[SR_MixLogNorm, histData, agediffsBinCounts, logSR_MixNorm, logSRbinCounts] = plotSRandResHistograms(nSRcounts, x, coreSubsetLogical, weightbydepthQ, weightRepDP, weightRepInflator, components, regularizationValue, subsetName, plotQ)
%%% This function takes some normalised sedimentation rate data in cell
%%% format, combines the counts into a single array, applies the weighting
%%% and then fits a mixture log normal to the result. It will also plot the
%%% histogram of the logarithm data with the mixture normal on top, show a
%%% QQ plot of the logarithm data, show the histogram of the data with the
%%% log mix normal on top, and plot a histogram of the resolution of the
%%% data

%% Combine all counts into one array
%set up arrays to be concatenated into
nSRcountsArray = ones(4,1);

if weightbydepthQ == 0
    %Concatenate all nSRcounts that are in the desired subset, choosing
    %only 1000 possible histories from each core (this allows the different
    %scenarios to be accounted for without letting them overweight the influence of that
    %core).
    for i = 1:length(nSRcounts)
        if coreSubsetLogical(i) == 1
            corenSRcounts = nSRcounts{i};
            numHistories = length(corenSRcounts);
            chosenHistoriesI = randi(numHistories, 1,1000);
            chosenHistories = corenSRcounts(:,chosenHistoriesI);
            nSRcountsArray = cat(2, nSRcountsArray, chosenHistories);
        end
    end
else
    %Concatenate all nSRcounts that are in the desired subset
    for i = 1:length(nSRcounts)
        if coreSubsetLogical(i) == 1
            nSRcountsArray = cat(2, nSRcountsArray, nSRcounts{i});
        end
    end
end

%Remove the ones that were used to set up arrays
nSRcountsArray = nSRcountsArray(:,2:end);

%% Fit Mix Log Norm
%Convert the weighted nSR counts to a single dimension array of counts
%where the number of counts is representative of their weighting
nSR = nSRcountsArray(1,:)'; %nSR data
depthWeights = nSRcountsArray(2,:);  %weightings
agedifferences = nSRcountsArray(4,:);
nSRclean = nSR(~isnan(nSR)); %Remove NaNs that separate cores and runs
dWeightsclean = depthWeights(~isnan(nSR)); %Remove NaNs that separate cores and runs
agedifferencesclean = agedifferences(~isnan(nSR));
if weightbydepthQ == 0
    data = nSRclean; 
else
    data = makeWeightedReplicates(nSRclean, dWeightsclean, weightRepDP, weightRepInflator);
end
dataLog = log(data);

%Calculate the MixLogNorm from these counts
[a,b] = size(data);
if a>b
    data = data';
end
[SR_MixLogNorm, logSR_MixNorm, gmfit] = fitMixLogNorm(data, x, components, regularizationValue);

%% count how many estimates of nSR
numbernSRcounts = length(nSRclean);

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
agediffsBinEdges = 0:500:10000;

%Calculate bin counts (using weighting for SR)
SRbinCounts = makeWeightedBinCounts(nSRcountsArray(1,:), nSRcountsArray(2,:), SRbinEdges);
histData = data;
logSRbinCounts = makeWeightedBinCounts(log(nSRclean), dWeightsclean, logSRbinEdges);
agediffsBinCounts = makeWeightedBinCounts(agedifferencesclean, dWeightsclean, agediffsBinEdges);

%% Compare histogram of "weighted data" to weighted histogram of data
figure;
subplot(1,2,1)
histogram("BinCounts", SRbinCounts, "BinEdges", SRbinEdges)
xlim([0 6])
title("Weighted by summing weights")
subplot(1,2,2)
histogram(data, "BinEdges", SRbinEdges)
xlim([0 6])
title("Weighted by replicating data")
% 
%% compare histogram of weighted log data to weighted histogram of log data
figure;
subplot(1,2,1)
histogram("BinCounts", logSRbinCounts, "BinEdges", logSRbinEdges)
xlim([-4 4])
title("Weighted by summing weights")
subplot(1,2,2)
histogram(dataLog, "BinEdges", logSRbinEdges)
xlim([-4 4])
title("Weighted by replicating data")
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
% 

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