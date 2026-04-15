function[countsArray] = countsCell2Array(countsCell, subsetLogical)
% countsCell2Array  Concatenate nSR matrices from a cell array into one array.
%
% Selects the cells of countsCell identified by subsetLogical and
% horizontally concatenates their contents into a single 3-row nSR matrix.
% Cells that are not selected (subsetLogical == 0) are skipped entirely.
%
% INPUTS
%   countsCell    - (cell array) One cell per core, each containing a
%                   3-row nSR matrix (see README "Internal Data Formats").
%                   NaN in the first row of a column marks a run header;
%                   data columns carry nSR, depth-difference, and
%                   age-difference values in rows 1–3 respectively.
%   subsetLogical - (logical vector) Same length as countsCell. A value of
%                   1 includes that cell's matrix in the output; 0 skips it.
%
% OUTPUT
%   countsArray   - (3 × M numeric matrix) All selected nSR matrices
%                   concatenated horizontally into a single array.
%
% See also: ARfitdists, makeWeightedBinCounts, makeWeightedReplicates

%Initialise array to be added to
countsArray = ones(3,1);

%Concatenate all nSRcounts that are in the desired subset
for i = 1:length(countsCell)
    if subsetLogical(i) == 1
        countsArray = cat(2, countsArray, countsCell{i});
    end
end
%Remove the ones that were used to set up arrays
countsArray = countsArray(:,2:end);