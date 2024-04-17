function[binCounts] = makeWeightedBinCounts(countsArray, weightingArray, binedges)
k = 0;
%For each bin...
binCounts = nan(1,length(binedges)-1);
for i = 1:length(binedges)-1
    %Find index of values that fit in the bin
    ix = (countsArray>= binedges(i) & countsArray < binedges(i+1));
    k = k+1;
    %Find counts for that bin by counting the weights (depth intervals)
    binCounts(k) = sum(weightingArray(ix));
    ix = []; %reset index vector
end
end