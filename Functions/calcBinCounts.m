function[obsCD, expCounts1, expCounts2] = calcBinCounts(data, binEdges, desiredSum, pdfV1,pdfV2)
%Find how many observations are in each bin
obsC  = histcounts(data, binEdges);
divisor = numel(data)/desiredSum;
obsCD = round(obsC./divisor);

%Get the expected counts given pdfVEC1
cdf1_at_edges = interp1(pdfV1.x, pdfV1.cdf_x, binEdges(2:end-1));
expProbs1 = diff([0, cdf1_at_edges, 1]);
expCounts1 = expProbs1.*desiredSum;

%Get the expected counts, given pdfVEC2
cdf2_at_edges = interp1(pdfV2.x, pdfV2.cdf_x, binEdges(2:end-1));
expProbs2 = diff([0, cdf2_at_edges, 1]);
expCounts2 = expProbs2.*desiredSum;
end