function[] = plotSRandResHistograms(nSRcounts, agediffs, num14cpairs, coreSubsetLogical, fignumber, colour, subsetName)

%% Combine all counts into one array
%set up arrays to be concatenated into
nSRcountsArray = ones(2,1);
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
%where the number of counts is representative of their weight
X = round(nSRcountsArray(1,:),2);
Y = nSRcountsArray(2,:);
[X_u, IA, IC] = unique(X);
Y_u = accumarray(IC,Y);
Y_uR = round(Y_u);
data = repelem(X_u,Y_uR);

%Calculate the MixLogNorm from these counts
nSRdataAsCounts = data;
SR_MixLogNorm = fitMixLogNorm(nSRdataAsCounts);

%% Calculate Histogram Bin Counts
%Define bins edges for histogram
SRbinEdges = 0:0.1:6;
agediffsBinEdges = 0:500:10000;

%Calculate bin counts (using weighting for SR)
SRbinCounts = makeWeightedBinCounts(nSRcountsArray(1,:), nSRcountsArray(2,:), SRbinEdges);
agediffsBinCounts = makeWeightedBinCounts(agediffsArray, ones(1,length(agediffsArray)), agediffsBinEdges);

%Plot histograms of SR and of agediffs
figure(fignumber)
subplot(1,2,1)
yyaxis right
plot(SR_MixLogNorm(:,1), SR_MixLogNorm(:,2), 'k-')
hold on
yyaxis left
histogram("BinCounts", SRbinCounts, "BinEdges", SRbinEdges);
xlim([0 6])
xlabel("Normalised Sed Rate")
ylabel("Counts")

subplot(1,2,2)
histogram("BinCounts", agediffsBinCounts, "BinEdges", agediffsBinEdges)
xlabel("Age Diff (yrs)")
ylabel("Counts")
title(subsetName + " with " + num2str(sum(subset14cpairs)) + " age pairs")
end