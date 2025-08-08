function[obsCD, expCountsArr] = calcBinCountsFH(data, binEdges, desiredSum, cdfFH)
%Find how many observations are in each bin
obsC  = histcounts(data, binEdges);
divisor = numel(data)/desiredSum;
obsCD = round(obsC./divisor);

%Get expected counts for each pdf
numcdfFH = 1;
expCountsArr = NaN(numcdfFH, length(binEdges)-1);
for i = 1:numcdfFH
    cdfi = cdfFH; %Choose pdf
    cdfAtEdges = cdfFH(binEdges(2:end-1)); %Get cdf values at each internal binEdge
    expProbsi = diff([0, cdfAtEdges, 1]); %Find difference of cdf values at each binEdge (bottom must be 0, top must be 1), this provides what fraction of probability is within this bin
    expCountsi = expProbsi.*desiredSum; % Multiply this by desired sum
    expCountsArr(i,:) = expCountsi; %Save results to array
end
end