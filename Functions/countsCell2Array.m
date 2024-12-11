function[countsArray] = countsCell2Array(countsCell, subsetLogical)

%Initialise array to be added to
countsArray = ones(4,1);

%Concatenate all nSRcounts that are in the desired subset
for i = 1:length(countsCell)
    if subsetLogical(i) == 1
        countsArray = cat(2, countsArray, countsCell{i});
    end
end
%Remove the ones that were used to set up arrays
countsArray = countsArray(:,2:end);