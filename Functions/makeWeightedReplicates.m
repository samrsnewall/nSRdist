function[repData] = makeWeightedReplicates(data, weight, dataRoundingDP, weightingInflationMultiplier)
%%% This function takes some data and some weighting and creates
%%% the variable repData, which is a single variable that approximates the
%%% weighting of the data by replicating each x by approximately it's
%%% weighting. For example, if x=3 has a weighting y=4, then the value 3
%%% will be repeated 4 times in repData.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(data) ~= length(weight)
    warning("Weights array not same size as data array - weights may not be properly aligned")
end

 data           = round(data,dataRoundingDP);                   %The data are rounded so that there are less unique values to help count weightings
 data           = data(data ~=0);                               %Zeros are removed (choice for use in gamma fitting);
 weight         = weight(data ~=0);
[data_u,~, IC]  = unique(data);                                 %The unique values of the data are found, with their indices
 weight_u       = accumarray(IC,weight);                        %The weighting of each unique value is combined
 weight_uR      = round(weight_u.*weightingInflationMultiplier);%The weighting of each unique value is converted to an integer
 repData        = repelem(data_u,weight_uR);                    %The unique values are replicated according to their integer weighting
end