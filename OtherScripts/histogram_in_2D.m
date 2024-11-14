%% Combine all counts into one array
%set up arrays to be concatenated into
nSRcountsArray = ones(4,1);
%agediffsArray = ones(1,1);
%subset14cpairs = zeros(1,length(num14cpairs));
%Concatenate all nSRcounts that are in the desired subset
for i = 1:length(nSRcounts)
    if coreSubsetLogical(i) == 1
    nSRcountsArray = cat(2, nSRcountsArray, nSRcounts{i});
    %agediffsArray = cat(2, agediffsArray, agediffs{i});
    %subset14cpairs(i) = num14cpairs(i);
    end
end
%Remove the ones that were used to set up arrays
nSRcountsArray = nSRcountsArray(:,2:end);
%agediffsArray = agediffsArray(:,2:end);

%% Remove all nans from the arrays
nSR = nSRcountsArray(1,:)'; %nSR data
depthWeights = nSRcountsArray(2,:);  %weightings
agedifferences = nSRcountsArray(4,:);
nSRclean = nSR(~isnan(nSR)); %Remove NaNs that separate cores and runs
dWeightsclean = depthWeights(~isnan(nSR)); %Remove NaNs that separate cores and runs
agedifferencesclean = agedifferences(~isnan(nSR));

lognSRclean = log(nSRclean);

%% Show a 2D histogram
[N, Xedges, Yedges] = histcounts2(nSRclean, agedifferencesclean');
imagesc(N);
title("nSR")

[N, Xedges, Yedges] = histcounts2(1./nSRclean, agedifferencesclean');
imagesc(N);
title("inverse nSR")

