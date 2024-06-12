function[] = plotSRandResHistograms(nSRcounts, agediffs, num14cpairs, coreSubsetLogical, ~, colour, subsetName)

%% Combine all counts into one array
%set up arrays to be concatenated into
nSRcountsArray = ones(3,1);
agediffsArray = ones(1,1);
subset14cpairs = zeros(1,length(num14cpairs));
%Concatenate all nSRcounts that are in the desired subset
for i = 1:length(nSRcounts)
    if coreSubsetLogical(i) == 1
    nSRcountsArray = cat(2, nSRcountsArray, nSRcounts{i});
    agediffsArray = cat(2, agediffsArray, agediffs{i});
    subset14cpairs(i) = num14cpairs(i);
    end
end
%Remove the ones that were used to set up arrays
nSRcountsArray = nSRcountsArray(:,2:end);
agediffsArray = agediffsArray(:,2:end);

%% Fit Mix Log Norm
%Convert the weighted nSR counts to a single dimension array of counts
%where the number of counts is representative of their weighting
X = round(nSRcountsArray(1,:),2)';
Y = nSRcountsArray(2,:);
[X_u, ~, IC] = unique(X);
Y_u = accumarray(IC,Y);
X_u = X_u(~isnan(X_u)); %Remove NaNs that separate cores and runs
Y_u = Y_u(~isnan(X_u)); %Remove NaNs that separate cores and runs
Y_uR = round(Y_u.*100); %%%%%%%%%% NEED TO DO A SENSITIVITY TEST TO THIS! IF I DON'T MULTIPY Y_u BY 100 THEN MUCH DATA HAS ITS WEIGHTING ROUNDED DOWN TO 0 WEIGHT... TRY 
data = repelem(X_u,Y_uR);
dataLog = log(data);

%Calculate the MixLogNorm from these counts
nSRdataAsCounts = data;
[SR_MixLogNorm, logSR_MixNorm, ~, ~] = fitMixLogNorm(nSRdataAsCounts, 2);

%% Calculate mean and variance of data
dataMeanLinear = geomean(nSRdataAsCounts);
dataVarLinear = var(nSRdataAsCounts);
dataMeanLog = mean(dataLog);
dataVarLog = var(dataLog);

%% Calculate Histogram Bin Counts
%Define bins edges for histogram
SRbinEdges = 0:0.1:6;
agediffsBinEdges = 0:500:10000;

%Calculate bin counts (using weighting for SR)
SRbinCounts = makeWeightedBinCounts(nSRcountsArray(1,:), nSRcountsArray(2,:), SRbinEdges);
agediffsBinCounts = makeWeightedBinCounts(agediffsArray, ones(1,length(agediffsArray)), agediffsBinEdges);



%% Plot to see how data compare to estimated distributions
figure
subplot(2,2,1)
yyaxis left
histogram(dataLog)
yyaxis right
plot(logSR_MixNorm(:,1), logSR_MixNorm(:,2))
xlabel("Log SR")
title("Data Mean = " + num2str(dataMeanLog) + "; Data Var = " + num2str(dataVarLog))
subplot(2,2,2)
qqplot(dataLog)
xlabel("Log SR Data")
subplot(2,2,3)
yyaxis left
histogram("BinCounts", SRbinCounts, "BinEdges", SRbinEdges)
yyaxis right
plot(SR_MixLogNorm(:,1), SR_MixLogNorm(:,2))
xlabel("SR")
xlim([0 6])
title("Data Mean = " + num2str(dataMeanLinear) + "; Data Var = " + num2str(dataVarLinear))
subplot(2,2,4)
histogram("BinCounts", agediffsBinCounts, "BinEdges", agediffsBinEdges)
xlabel("Age Diff (yrs)")
ylabel("Counts")
title(num2str(sum(subset14cpairs)) + " age pairs")
sgtitle(subsetName)

figure
subplot(1,2,1)
yyaxis left
histogram("BinCounts", SRbinCounts, "BinEdges", SRbinEdges, 'FaceColor', colour, 'FaceAlpha', 0.5)
yyaxis right
plot(SR_MixLogNorm(:,1), SR_MixLogNorm(:,2), 'Color', colour, 'LineWidth', 2)
xlabel("SR")
xlim([0 6])
title("Data Mean = " + num2str(dataMeanLinear) + "; Data Var = " + num2str(dataVarLinear))
subplot(1,2,2)
histogram("BinCounts", agediffsBinCounts, "BinEdges", agediffsBinEdges,'FaceColor', colour, 'FaceAlpha', 0.5)
xlabel("Age Diff (yrs)")
ylabel("Counts")
title(num2str(sum(subset14cpairs)) + " age pairs")
sgtitle(subsetName)

end