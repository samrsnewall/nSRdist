function[repData] = makeWeightedReplicates(data, weight, dataRoundingDP, weightsMultiplier)
%%% This function takes some data and some weighting and creates
%%% the variable repData, which is a single variable that approximates the
%%% weighting of the data by replicating each x by approximately it's
%%% weighting. For example, if x=3 has a weighting y=4, then the value 3
%%% will be repeated 4 times in repData.

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