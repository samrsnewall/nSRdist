function[binCounts] = makeWeightedBinCounts(countsArray, weightingArray, binedges)
% makeWeightedBinCounts  Compute weighted histogram bin counts.
%
% For each bin defined by binedges, sums the weights (from weightingArray)
% of all data values (from countsArray) that fall within that bin. This
% produces a weighted histogram where each observation contributes its
% weight rather than a count of 1. Values exactly equal to a lower edge
% are included; values equal to an upper edge fall into the next bin
% (half-open interval [lower, upper)).
%
% INPUTS
%   countsArray    - (numeric vector) Data values to be binned
%   weightingArray - (numeric vector) Weight for each element of
%                    countsArray; must be the same length as countsArray
%   binedges       - (numeric vector, length B+1) Edges defining B bins
%
% OUTPUT
%   binCounts      - (1 × B numeric vector) Weighted count for each bin
%
% See also: ARfitdists, IRfitdists, countsCell2Array

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