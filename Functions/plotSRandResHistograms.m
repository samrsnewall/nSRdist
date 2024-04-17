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

%% Plot histograms
%Define bins edges for histogram
SRbinEdges = 0:0.01:6;
agediffsBinEdges = 0:500:10000;

%Calculate bin counts (using weighting for SR)
SRbinCounts = makeWeightedBinCounts(nSRcountsArray(1,:), nSRcountsArray(2,:), SRbinEdges);
agediffsBinCounts = makeWeightedBinCounts(agediffsArray, ones(1,length(agediffsArray)), agediffsBinEdges);

%Plot histograms of SR and of agediffs
figure(fignumber)
subplot(1,2,1)
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