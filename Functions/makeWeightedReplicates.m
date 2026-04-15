function[repData] = makeWeightedReplicates(data, weight, dataRoundingDP, weightsMultiplier)
% makeWeightedReplicates  Approximate a continuous weighting by replicating
%                         data values according to their weights.
%
% Distribution-fitting functions (e.g. fitgmdist, gamfit) require a vector
% of individual observations and cannot accept histogram bin counts or
% continuous weights directly. This function bridges that gap by converting
% a weighted dataset into an unweighted replicate dataset that approximates
% the same weighting.
%
% The procedure is:
%   1. Round data to dataRoundingDP decimal places, reducing the number of
%      unique values.
%   2. Discard any values that round to exactly zero (data are expected to
%      be strictly positive).
%   3. For each unique rounded value, accumulate the total weight across all
%      observations that share that value.
%   4. Multiply each accumulated weight by weightsMultiplier and round to
%      the nearest integer, giving an integer replication count.
%   5. Use repelem to replicate each unique value by its integer count.
%
% The resulting repData vector has total length equal to the sum of all
% integer replication counts. Values with very small accumulated weights
% may round to zero and be silently dropped; increase weightsMultiplier to
% reduce this effect if necessary.
%
% INPUTS
%   data             - (numeric vector) Data values to be replicated
%                      (strictly positive; values rounding to 0 are removed)
%   weight           - (numeric vector) Non-negative weight for each element
%                      of data; must be the same length as data
%   dataRoundingDP   - (integer) Number of decimal places to which data
%                      values are rounded before accumulation
%   weightsMultiplier - (scalar) Factor applied to accumulated weights
%                       before rounding to integers; increase to reduce
%                       precision loss from rounding small weights to zero
%
% OUTPUT
%   repData - (numeric vector) Replicated data approximating the weighted
%             distribution; length equals sum of all integer replication counts
%
% See also: ARfitdists, IRfitdists, makeWeightedBinCounts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(data) ~= length(weight)
    warning("Weights array not same size as data array - weights may not be properly aligned")
end

 data           = round(data,dataRoundingDP);                   %The data are rounded (to some decimal places) so that there are less unique values to help count weightings
 data           = data(data ~=0);                               %Remove any values that rounded down to zero (data should be strictly positive);
 weight         = weight(data ~=0);
[data_u,~, IC]  = unique(data);                                 %The unique values of the data are found, with their indices
 weight_u       = accumarray(IC,weight);                        %The weighting of each unique value is combined
 weight_uR      = round(weight_u.*weightsMultiplier);           %The weighting of each unique value is multiplied by some value (useful to avoid rounding to 0) and then rounded to an integer
 repData        = repelem(data_u,weight_uR);                    %The unique values are replicated according to their integer weighting



 % %Test to see if the approach above makes much of a difference to the very
 % %simple approach of
 % weight_R = round(weight.*weightsMultiplier);
 % repData2 = repelem(data, weight_R);
 % 
 % binEdges = 0:0.05:10;
 % 
 % figure;
 % subplot(2,1,1)
 % histogram(repData, 'BinEdges', binEdges)
 % title("repData")
 % subplot(2,1,2)
 % histogram(repData2, 'BinEdges', binEdges)
 % title("repData2")
 % 
 % if any(weight_uR == 0) || any(weight_R == 0)
 %     lostWeight1 = sum(weight_u(weight_uR == 0));
 %     lostWeight2 = sum(weight(weight_R == 0));
 %     warning("repData lost " + num2str(lostWeight1) + " weight ; repData2 lost " + num2str(lostWeight2) + " weight")
 % end

end